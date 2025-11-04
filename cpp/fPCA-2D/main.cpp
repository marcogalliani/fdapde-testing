//
// TEMPLATE: Minimal R ↔ C++ JSON interface
// Reads paths and options, loads data, calls `fit_model`, saves output
//

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

// Solver router
using fpca_solver_variant = std::variant<
  fpca_power_solver,
  fpca_subspace_solver,
  fpca_subspace_experimental_solver,
  fpca_direct_solver
>;
fpca_solver_variant get_fpca_solver(const std::string& solver_name) {
  if (solver_name == "subspace") return fpca_subspace_solver();
  else if (solver_name == "subspace_fpc_spec") return fpca_subspace_experimental_solver();
  else if (solver_name == "sequential") return fpca_power_solver();
  else if (solver_name == "direct") return fpca_direct_solver();
  else throw std::invalid_argument("Unknown solver: " + solver_name);
}

// Example function
auto fit_model(Triangulation<2,2> D,
                const matrix_t& X,
                const matrix_t& locs,
                const std::vector<double>& lambda_grid,
                const std::string& solver_name, 
                const int n_comp,
                bool compute_mean = false) {
  
  std::cout << "Fit: " << solver_name << std::endl;
  std::cout << "- Running solver: " << solver_name << std::endl;
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
  
  // Initialize the model
  fPCA model("X", data, fe_ls_elliptic(a, F));
  
  int mean = 0x0;
  if(compute_mean) mean = ComputeMean;

  // Select the fPCA solver according to solver_name
  std::visit(
    [&](auto&& solver){
      model.fit(n_comp, lambda_grid, ComputeRandSVD | OptimizeGCV | mean, solver);
    },
    get_fpca_solver(solver_name)
  );
  
  return model;
}

// Main
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
  std::string solver_name = jroot["options"].value("solver", "default_solver");
  std::vector<double> lambda_grid = jroot["options"].at("lambda_grid").get<std::vector<double>>();
  double n_comp = jroot["options"].value("n_comp", 3);
  bool mean = jroot["options"].value("mean", false);
  
  std::cout << "Options:" << std::endl;
  std::cout << "- Solver name: " << solver_name << std::endl;
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
  auto model = fit_model(D, X, locs, lambda_grid, solver_name, n_comp, mean);
  
  // Post-processing
  matrix_t rec_X = model.S()*model.F().transpose();
  // rec_X = rec_X.rowwise() + model.center().transpose();
  matrix_t rec_X_locs = model.S()*model.Fn().transpose();
  // rec_X_locs = rec_X_locs.rowwise() + model.center_locs().transpose();
  
  // Save results ----
  if(mean){
    write_csv(path_results + "center.csv", model.center());
    write_csv(path_results + "center_locs.csv", model.center_locs());
  }
  write_csv(path_results + "loadings.csv", model.F());
  write_csv(path_results + "loadings_locs.csv", model.Fn());
  write_csv(path_results + "scores.csv", model.S());
  write_csv(path_results + "reconstruction.csv", rec_X);
  write_csv(path_results + "reconstruction_at_locs.csv", rec_X_locs);
  write_csv(path_results + "lambda.csv", model.lambda());
  write_csv(path_results + "gcv_scores.csv", model.gcv_scores());
  //write_csv(path_results + "var_pct.csv", model.var_explained());
  //write_csv(path_results + "smoothed_data.csv", model.smoothed_data());
  
  // Ensure the file ends with a newline (optional cleanup)
  std::ofstream fix_newline(path_results + "lambda.csv", std::ios::app);
  fix_newline << std::endl;
  
  std::cout << "Results:" << std::endl;
  std::cout << "- Results written to: " << path_results << "output.csv\n";
  std::cout << std::endl;
  
  return 0;
}