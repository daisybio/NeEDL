//
// Created by juli on 10.05.23.
//

#ifndef GENEPISEEKER_SHINYAPPLAUNCHER_HPP
#define GENEPISEEKER_SHINYAPPLAUNCHER_HPP

#include "Job.hpp"
#include "InstanceLoader.hpp"

namespace epi {

    class ShinyAppLauncher : public Job {
    public:
        ShinyAppLauncher(std::shared_ptr<InstanceLoader> instance_loader);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        std::shared_ptr<InstanceLoader> instance_loader;
    };

} // epi

#ifdef HEADER_ONLY
#include "ShinyAppLauncher.cpp"
#endif

#endif //GENEPISEEKER_SHINYAPPLAUNCHER_HPP
