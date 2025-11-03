#include <fdaPDE/models.h>
using namespace fdapde;

#include "../include/json.hpp"
using nlohmann::json;

#include <Eigen/Dense>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <variant>

// Type aliases
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using sparse_matrix_t = Eigen::SparseMatrix<double>;

struct sMVPCA{
    matrix_t S_, F_, Fn_;
    vector_t lambda_;
    // observers
    const matrix_t& S() { return S_;}
    const matrix_t& F() { return F_;}
    const matrix_t& Fn() { return Fn_;}
    const vector_t& lambda() { return lambda_;} 
};

// Example function
auto fit_model(Triangulation<2,2> D,
                const matrix_t& X,
                const matrix_t& locs,
                const std::vector<double>& lambda_grid,
                const int n_comp) {
    std::cout << "Fit presmoothing-PCA " << std::endl;
    std::cout << "- Lambda grid has " << lambda_grid.size() << " values.\n";
    std::cout << std::endl;
    
    // Physics (isotropic Laplacian)
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    ZeroField<2> u;
    auto a = integral(D)(dot(grad(f), grad(v)));
    auto F = integral(D)(u * v);
    
    // GeoFrame
    GeoFrame data(D);
    auto& l = data.insert_scalar_layer<POINT>("locs_layer", locs);
    l.load_blk("X", X.transpose());
  
    sMVPCA model;
    // (1) pre-smoothing step
    // modeling
    FunctSmoother m("X", data, fe_ls_elliptic(a, F));
    // fitting
    m.fit(lambda_grid, OptimizeGCV);
    model.lambda_ = m.lambda();
    // (2) Generalized PCA
    // -> Cholesky of mass matrix
    Eigen::SimplicialLLT<sparse_matrix_t> llt(m.mass());
    // -> SVD of smoothed data left-multiplied by the cholesky factor of the mass matrix
    Eigen::JacobiSVD<matrix_t> svd;
    sparse_matrix_t cholesky_factor = llt.matrixL();
    cholesky_factor = llt.permutationPinv() * cholesky_factor;
    svd.compute(m.smoothed_data() * cholesky_factor, Eigen::ComputeThinU | Eigen::ComputeThinV);

    int n_dof = m.smoothed_data().cols();
    matrix_t Id = matrix_t::Identity(n_dof, n_dof);
    model.S_ = svd.matrixU().leftCols(n_comp)*svd.singularValues().head(n_comp).asDiagonal();
    model.F_ = llt.permutationPinv()*llt.matrixL().solve(Id).transpose()*svd.matrixV().leftCols(n_comp);
    model.Fn_ = m.Psi() * model.F_;

    return model;
}

int main(int argc, char* argv[]) {
    std::cout << std::endl;
  
    // Check for argument
    if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <params.json>" << std::endl;
    return 1;
    }

    std::string params_path = argv[1];
    std::cout << "Reading parameters from: " << params_path << std::endl;

    // Load JSON
    std::ifstream input(params_path);
    if (!input) {
    throw std::runtime_error("Cannot open params.json in current directory.");
    }
    json jroot = json::parse(input);

    // Close and delete the file
    input.close();
    std::filesystem::remove(params_path);

    // Extract paths
    std::string path_mesh = "../../" + jroot["path_list"].value("mesh",    "./mesh/");
    std::string path_data = "../../" + jroot["path_list"].value("data",    "./data/");
    std::string path_results = "../../" + jroot["path_list"].value("results", "./results/");

    std::cout << std::endl;
    std::cout << "Paths:" << std::endl;
    std::cout << "- Mesh: " << path_mesh << std::endl;
    std::cout << "- Data: " << path_data << std::endl;
    std::cout << "- Results: " << path_results << std::endl;
    std::cout << std::endl;

    // Extract options
    std::vector<double> lambda_grid = jroot["options"].at("lambda_grid").get<std::vector<double>>();
    double n_comp = jroot["options"].value("n_comp", 3);

    std::cout << "Options:" << std::endl;
    std::cout << "- Lambda grid: [" << lambda_grid[0] << " " << lambda_grid[1] << " ... " << lambda_grid.back() << "]" << std::endl;
    std::cout << "- N. comp: " << n_comp << std::endl;
    std::cout << std::endl;

    // Load geometry
    Triangulation<2,2> D(
        path_mesh + "points.csv",
        path_mesh + "elements.csv",
        path_mesh + "boundary.csv",
        true, true
    );

    // Load data
    matrix_t X = read_csv<double>(path_data + "X.csv").as_matrix();
    matrix_t locs = read_csv<double>(path_data + "locs.csv").as_matrix();

    std::cout << "Loaded data:" << std::endl;
    std::cout << "- X(" << X.rows() << ", " << X.cols() << ")" << std::endl;
    std::cout << "- locs(" << locs.rows() << ", " << locs.cols() << ")" << std::endl;
    std::cout << std::endl;

    // Fit the model
    auto model = fit_model(D, X, locs, lambda_grid, n_comp);
    
    // Post-processing
    // rec_X = rec_X.rowwise() + model.center().transpose();
    matrix_t rec_X_locs = model.S()*model.Fn().transpose();
    // rec_X_locs = rec_X_locs.rowwise() + model.center_locs().transpose();
    
    // Save results ----
    write_csv(path_results + "loadings_locs.csv", model.Fn());
    write_csv(path_results + "loadings.csv", model.F());
    write_csv(path_results + "scores.csv", model.S());
    write_csv(path_results + "reconstruction_at_locs.csv", rec_X_locs);
    write_csv(path_results + "lambda.csv", model.lambda());
    
    // Ensure the file ends with a newline (optional cleanup)
    std::ofstream fix_newline(path_results + "lambda.csv", std::ios::app);
    fix_newline << std::endl;
    
    std::cout << "Results:" << std::endl;
    std::cout << "- Results written to: " << path_results << "output.csv\n";
    std::cout << std::endl;

    return 0;
}
