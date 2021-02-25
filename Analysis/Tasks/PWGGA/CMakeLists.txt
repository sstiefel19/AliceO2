o2_add_dpl_workflow(ga-gammajet
                    SOURCES gammaJet.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisJets O2::AnalysisDataModel O2::DetectorsBase O2::AnalysisCore
                    COMPONENT_NAME Analysis)
