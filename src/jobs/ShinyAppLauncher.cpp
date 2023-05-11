//
// Created by juli on 10.05.23.
//

#include "ShinyAppLauncher.hpp"

namespace epi {
    void ShinyAppLauncher::run(std::shared_ptr<DataModel> data) {
        std::string dataset_path = instance_loader->getInputFilePath();
        std::string results_path = data->outputDirectory->get_output_directory();

        std::string dataset_path_relative = boost::filesystem::relative(dataset_path, results_path).string();
        std::string dataset_path_absolute = boost::filesystem::absolute(dataset_path).string();

        std::string launcher_linux_path = data->outputDirectory->get_output_directory() + "SHOW_RESULTS_LINUX.sh";
        std::ofstream launcher_linux (launcher_linux_path);
        launcher_linux << "#!/bin/bash\n"
                          "if [ -f \"" << dataset_path_relative << "\" ]; then\n"
                          "    dataset=\"" << dataset_path_relative << "\"\n"
                          "else\n"
                          "    dataset=\"" << dataset_path_absolute << "\"\n"
                          "fi\n"
                          "dataset=$(realpath \"$dataset\")\n"
                          "results=$(realpath ./)\n"
                          "dbsnp_internal=/NeEDL/data/dbSNP/inc_pseudogenes/snps_restruc_full_inc_pseudo.csv\n"
                          "biogrid_internal=/NeEDL/data/BIOGRID/BIOGRID-ORGANISM-Homo_sapiens-current.tab2.txt\n"
                          "dataset_internal=/app/data/dataset.json\n"
                          "results_internal=/app/data/results/\n"
                          "user_id=$(id -u)\n"
                          "group_id=$(id -g)\n"
                          "docker run --rm -v \"$dataset:$dataset_internal:Z\" \\\n"
                          "    -v \"$results:$results_internal:Z\" \\\n"
                          "    --user $user_id:$group_id \\\n"
                          "    --workdir /app \\\n"
                          "    --network host \\\n"
                          "    ghcr.io/biomedbigdata/genepiseeker_dev:master Rscript app.R \\\n"
                          "    --results $results_internal \\\n"
                          "    --dataset $dataset_internal \\\n"
                          "    --dbSNP $dbsnp_internal \\\n"
                          "    --biogrid $biogrid_internal\n";

        launcher_linux.close();

        boost::filesystem::permissions(launcher_linux_path, boost::filesystem::perms::add_perms | boost::filesystem::perms::owner_exe);

    }

    ShinyAppLauncher::ShinyAppLauncher(std::shared_ptr<InstanceLoader> instance_loader) {
        this->instance_loader = instance_loader;
    }
} // epi