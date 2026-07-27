#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TAxis.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

struct ParameterResult {
    double value = 0.0;
    double error = 0.0;
};

struct BootstrapResult {
    int seed = -1;
    int status = -999;
    int cov_matrix_status = -999;
    double loglikelihood_min_value = 0.0;
    std::map<std::string, ParameterResult> parameters;
    std::string source_file;
};

std::string Trim(const std::string& text)
{
    const std::size_t first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";
    const std::size_t last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

bool ParseSeed(const std::string& line, int& seed)
{
    static const std::regex pattern(
        R"(^\s*seed\s*=\s*(-?\d+)\s*$)"
    );
    std::smatch match;

    if (!std::regex_match(line, match, pattern))
        return false;
    seed = std::stoi(match[1].str());
    return true;
}

bool ParseStatusLine(
    const std::string& line,
    int& status,
    int& cov_matrix_status,
    double& loglikelihood_min_value)
{
    static const std::regex pattern(
        R"(^\s*status\s*=\s*(-?\d+)\s+CovMatrixStatus\s*=\s*(-?\d+)\s+loglikehood_min_value\s*=\s*([-+0-9.eE]+)\s*$)"
    );

    std::smatch match;

    if (!std::regex_match(line, match, pattern))
        return false;

    status = std::stoi(match[1].str());
    cov_matrix_status = std::stoi(match[2].str());
    loglikelihood_min_value = std::stod(match[3].str());

    return true;
}

bool ParseParameterLine(
    const std::string& line,
    std::string& parameter_name,
    double& value,
    double& error)
{
    static const std::regex pattern(
        R"(^\s*([^:]+?)\s*:\s*([-+0-9.eE]+)\s*\+\-\s*([-+0-9.eE]+)\s*$)"
    );

    std::smatch match;

    if (!std::regex_match(line, match, pattern))
        return false;

    parameter_name = Trim(match[1].str());
    value = std::stod(match[2].str());
    error = std::stod(match[3].str());

    return true;
}

bool ReadBootstrapFile(
    const std::string& filename,
    std::string& reference_header,
    bool& reference_header_initialized,
    std::vector<BootstrapResult>& results,
    std::vector<std::string>& reference_parameter_names,
    bool& reference_parameters_initialized)
{
    std::ifstream input(filename);

    if (!input.is_open()) {
        std::cerr << "ERROR: cannot open " << filename << '\n';
        return false;
    }

    std::string line;

    if (!std::getline(input, line)) {
        std::cerr << "ERROR: empty file " << filename << '\n';
        return false;
    }

    const std::string current_header = line;

    if (!reference_header_initialized) {
        reference_header = current_header;
        reference_header_initialized = true;
    } else if (current_header != reference_header) {
        std::cerr << "ERROR: first line differs in file:\n"
                  << "  " << filename << '\n'
                  << "Reference:\n"
                  << "  " << reference_header << '\n'
                  << "Current:\n"
                  << "  " << current_header << '\n';

        return false;
    }

    BootstrapResult current_result;
    bool inside_block = false;
    bool status_found = false;

    auto close_current_block = [&]() -> bool {
        if (!inside_block)
            return true;

        if (!status_found) {
            std::cerr << "WARNING: missing status line for seed "
                      << current_result.seed << " in " << filename << '\n';
        }

        std::vector<std::string> current_parameter_names;

        for (const auto& [name, result] : current_result.parameters)
            current_parameter_names.push_back(name);

        if (!reference_parameters_initialized) {
            reference_parameter_names = current_parameter_names;
            reference_parameters_initialized = true;
        } else if (current_parameter_names != reference_parameter_names) {
            std::cerr << "ERROR: variable list differs for seed "
                      << current_result.seed << " in " << filename << '\n';

            std::cerr << "Expected variables:\n";

            for (const auto& name : reference_parameter_names)
                std::cerr << "  " << name << '\n';

            std::cerr << "Found variables:\n";

            for (const auto& name : current_parameter_names)
                std::cerr << "  " << name << '\n';

            return false;
        }

        results.push_back(current_result);

        return true;
    };

    while (std::getline(input, line)) {
        line = Trim(line);

        if (line.empty())
            continue;

        int seed = -1;

        if (ParseSeed(line, seed)) {
            if (!close_current_block())
                return false;

            current_result = BootstrapResult();
            current_result.seed = seed;
            current_result.source_file = filename;

            inside_block = true;
            status_found = false;

            continue;
        }

        if (!inside_block) {
            std::cerr << "WARNING: ignoring line outside a block in "
                      << filename << ":\n  " << line << '\n';
            continue;
        }

        int status = 0;
        int cov_matrix_status = 0;
        double loglikelihood = 0.0;

        if (ParseStatusLine(
                line,
                status,
                cov_matrix_status,
                loglikelihood)) {
            current_result.status = status;
            current_result.cov_matrix_status = cov_matrix_status;
            current_result.loglikelihood_min_value = loglikelihood;
            status_found = true;

            continue;
        }

        std::string parameter_name;
        double value = 0.0;
        double error = 0.0;

        if (ParseParameterLine(line, parameter_name, value, error)) {
            if (current_result.parameters.find(parameter_name) != current_result.parameters.end()) {
                std::cerr << "ERROR: duplicated parameter "
                          << parameter_name
                          << " for seed " << current_result.seed
                          << " in " << filename << '\n';

                return false;
            }

            current_result.parameters[parameter_name] = {value, error};

            continue;
        }

        std::cerr << "WARNING: unrecognized line in "
                  << filename << ":\n  " << line << '\n';
    }

    if (!close_current_block())
        return false;

    return true;
}

void PlotBootstrapResults(
    const char* directory = ".",
    const char* file_prefix = "bootstrapout",
    const char* output_prefix = "bootstrap_summary")
{
    namespace fs = std::filesystem;

    const std::string directory_name = directory;
    const std::string prefix = file_prefix;
    const std::string output_name = output_prefix;

    if (!fs::exists(directory_name)) {
        std::cerr << "ERROR: directory does not exist: "
                  << directory_name << '\n';
        return;
    }

    const std::regex filename_pattern(
        "^" + prefix + R"(_([0-9]+)\.txt$)"
    );

    std::vector<std::pair<int, std::string>> input_files;

    for (const auto& entry : fs::directory_iterator(directory_name)) {
        if (!entry.is_regular_file())
            continue;

        const std::string basename = entry.path().filename().string();
        std::smatch match;

        if (!std::regex_match(basename, match, filename_pattern))
            continue;

        const int file_number = std::stoi(match[1].str());

        input_files.push_back({
            file_number,
            entry.path().string()
        });
    }

    if (input_files.empty()) {
        std::cerr << "ERROR: no files matching "
                  << prefix << "_numero.txt found in "
                  << directory_name << '\n';
        return;
    }

    std::sort(
        input_files.begin(),
        input_files.end(),
        [](const auto& first, const auto& second) {
            return first.first < second.first;
        }
    );

    std::cout << "Found " << input_files.size()
              << " input files:\n";

    for (const auto& [number, filename] : input_files)
        std::cout << "  " << filename << '\n';

    std::vector<BootstrapResult> results;
    std::string reference_header;
    bool reference_header_initialized = false;

    std::vector<std::string> parameter_names;
    bool reference_parameters_initialized = false;

    for (const auto& [number, filename] : input_files) {
        if (!ReadBootstrapFile(
                filename,
                reference_header,
                reference_header_initialized,
                results,
                parameter_names,
                reference_parameters_initialized)) {
            std::cerr << "ERROR while reading " << filename << '\n';
            return;
        }
    }

    if (results.empty()) {
        std::cerr << "ERROR: no bootstrap blocks found.\n";
        return;
    }

    std::sort(
        results.begin(),
        results.end(),
        [](const BootstrapResult& first,
           const BootstrapResult& second) {
            return first.seed < second.seed;
        }
    );

    std::set<int> seeds;

    for (const BootstrapResult& result : results) {
        if (!seeds.insert(result.seed).second) {
            std::cerr << "ERROR: duplicated seed "
                      << result.seed
                      << " found in file "
                      << result.source_file << '\n';
            return;
        }
    }

    std::cout << "\nHeader check passed.\n";
    std::cout << "Total blocks: " << results.size() << '\n';
    std::cout << "Number of variables: "
              << parameter_names.size() << '\n';

    for (const std::string& name : parameter_names)
        std::cout << "  " << name << '\n';

    gStyle->SetOptStat(0);

    const std::string root_filename =
        output_name + ".root";

    const std::string pdf_filename =
        output_name + ".pdf";

    TFile output_file(root_filename.c_str(), "RECREATE");

    if (output_file.IsZombie()) {
        std::cerr << "ERROR: cannot create "
                  << root_filename << '\n';
        return;
    }

    TCanvas canvas("canvas", "Bootstrap results", 1200, 800);

    canvas.Print((pdf_filename + "[").c_str());

    for (const std::string& parameter_name : parameter_names) {
        TGraphErrors valid_graph;
        TGraphErrors invalid_graph;

        valid_graph.SetName(
            ("graph_" + parameter_name + "_valid").c_str()
        );

        invalid_graph.SetName(
            ("graph_" + parameter_name + "_invalid").c_str()
        );

        int valid_point = 0;
        int invalid_point = 0;

        for (const BootstrapResult& result : results) {
            const ParameterResult& parameter =
                result.parameters.at(parameter_name);

            const bool valid_fit =
                result.status == 0 &&
                result.cov_matrix_status == 3;

            if (valid_fit) {
                valid_graph.SetPoint(
                    valid_point,
                    result.seed,
                    parameter.value
                );

                valid_graph.SetPointError(
                    valid_point,
                    0.0,
                    parameter.error
                );

                ++valid_point;
            } else {
                invalid_graph.SetPoint(
                    invalid_point,
                    result.seed,
                    parameter.value
                );

                invalid_graph.SetPointError(
                    invalid_point,
                    0.0,
                    parameter.error
                );

                ++invalid_point;
            }
        }

        canvas.Clear();
        canvas.SetGrid();

        valid_graph.SetTitle(
            Form(
                "%s;Seed;%s",
                parameter_name.c_str(),
                parameter_name.c_str()
            )
        );

        valid_graph.SetMarkerStyle(20);
        valid_graph.SetMarkerSize(0.9);
        valid_graph.SetMarkerColor(kBlue + 1);
        valid_graph.SetLineColor(kBlue + 1);

        invalid_graph.SetMarkerStyle(24);
        invalid_graph.SetMarkerSize(1.0);
        invalid_graph.SetMarkerColor(kRed + 1);
        invalid_graph.SetLineColor(kRed + 1);

        if (valid_graph.GetN() > 0) {
            valid_graph.Draw("AP");
            invalid_graph.Draw("P SAME");
        } else if (invalid_graph.GetN() > 0) {
            invalid_graph.SetTitle(
                Form(
                    "%s;Seed;%s",
                    parameter_name.c_str(),
                    parameter_name.c_str()
                )
            );

            invalid_graph.Draw("AP");
        }

        TLegend legend(0.72, 0.78, 0.9, 0.9);

        if (valid_graph.GetN() > 0)
            legend.AddEntry(
                &valid_graph,
                "status = 0, CovMatrixStatus = 3",
                "lep"
            );

        if (invalid_graph.GetN() > 0)
            legend.AddEntry(
                &invalid_graph,
                "Other status",
                "lep"
            );

        legend.Draw();

        canvas.Print(pdf_filename.c_str());

        output_file.cd();
        valid_graph.Write();
        invalid_graph.Write();
    }

    TGraphErrors loglikelihood_valid;
    TGraphErrors loglikelihood_invalid;

    loglikelihood_valid.SetName("loglikelihood_valid");
    loglikelihood_invalid.SetName("loglikelihood_invalid");

    int valid_point = 0;
    int invalid_point = 0;

    for (const BootstrapResult& result : results) {
        const bool valid_fit =
            result.status == 0 &&
            result.cov_matrix_status == 3;

        if (valid_fit) {
            loglikelihood_valid.SetPoint(
                valid_point,
                result.seed,
                result.loglikelihood_min_value
            );

            ++valid_point;
        } else {
            loglikelihood_invalid.SetPoint(
                invalid_point,
                result.seed,
                result.loglikelihood_min_value
            );

            ++invalid_point;
        }
    }

    canvas.Clear();
    canvas.SetGrid();

    loglikelihood_valid.SetTitle(
        "Minimum likelihood value;Seed;loglikehood_min_value"
    );

    loglikelihood_valid.SetMarkerStyle(20);
    loglikelihood_valid.SetMarkerColor(kBlue + 1);
    loglikelihood_valid.SetLineColor(kBlue + 1);

    loglikelihood_invalid.SetMarkerStyle(24);
    loglikelihood_invalid.SetMarkerColor(kRed + 1);
    loglikelihood_invalid.SetLineColor(kRed + 1);

    if (loglikelihood_valid.GetN() > 0) {
        loglikelihood_valid.Draw("AP");
        loglikelihood_invalid.Draw("P SAME");
    } else {
        loglikelihood_invalid.SetTitle(
            "Minimum likelihood value;Seed;loglikehood_min_value"
        );

        loglikelihood_invalid.Draw("AP");
    }

    TLegend likelihood_legend(0.72, 0.78, 0.9, 0.9);

    if (loglikelihood_valid.GetN() > 0)
        likelihood_legend.AddEntry(
            &loglikelihood_valid,
            "Valid fit",
            "p"
        );

    if (loglikelihood_invalid.GetN() > 0)
        likelihood_legend.AddEntry(
            &loglikelihood_invalid,
            "Other status",
            "p"
        );

    likelihood_legend.Draw();

    canvas.Print(pdf_filename.c_str());
    canvas.Print((pdf_filename + "]").c_str());

    output_file.cd();
    loglikelihood_valid.Write();
    loglikelihood_invalid.Write();
    output_file.Close();

    std::cout << "\nCreated:\n"
              << "  " << pdf_filename << '\n'
              << "  " << root_filename << '\n';
}