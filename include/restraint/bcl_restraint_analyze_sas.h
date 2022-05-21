// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_RESTRAINT_ANALYZE_SAS_H_
#define BCL_RESTRAINT_ANALYZE_SAS_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_sas_scattering_data.h"
#include "bcl_restraint_sas_transformation.h"
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_model_data.h"
#include "score/bcl_score_sas_type.h"
#include "util/bcl_util_object_data_label.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeSas
    //! @brief performs analysis of protein or bcl model with experimental SAXS data
    //! @details can optionally display desired parameters
    //!
    //! @see @link example_restraint_analyze_saxs.cpp @endlink
    //! @author putnamdk
    //! @date Sep 5, 2014
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzeSas :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the postfix that will be appended to the filename that will hold this analysis
      std::string m_OutFilePostFix;

      //! the experimental SAXS data to compare computed profile to
      std::string m_ExperimentalDataFileName;

      //! the solvent accessible surface area file name to compute the hydration shell
      std::string m_SasaDataFileName;

      //! the computed saxs curve with experimental data file name
      std::string m_ComputedDataFileName;

      //! The c1 tuning parameter for excluded particle volume
      float m_ExcludedVolumeParameter;

      //! The c2 tuning parameter for hydration shell thickness
      float m_HydrationShellThickness;

      //! The maximum dimension of the p(r) from the experimental data set
      float m_MaximumDimension;

      //! true if search is to be performed for optimal c1 and c2 parameters
      bool m_OptimizeHydrationShellParameters;

      //! c1 and c2 search grid parameters
      float m_C1Min;
      float m_C1Max;
      float m_C2Min;
      float m_C2Max;
      float m_C1StepSize;
      float m_C2StepSize;

      //! the name of the scoring function to use to compare two SAXS curves
      score::SasType::ScoreFunctionEnum m_ScoreType;

      //! varible to control the state the data is presented
      //! i.e. scaled, normalized, log10, derivative
      SasTransformation m_Transform;

      //! varible to set default search grid
      //! m_C1Min( 0.80), m_C1Max( 1.20), m_C2Min( 0.0), m_C2Max( 4.0), m_C1StepSize( 0.005),  m_C2StepSize( 0.1)
      bool m_DefaultSearchGrid;

      //! switch to use cpu instead of gpu processing
      bool m_Cpu;

      //! flag to approximate Side Chains
      bool m_ShouldApproximateSideChains;

      bool m_ShouldApproximateLoops;

      //! flag to use SANS
      bool m_ShouldUseSans;

      //! value between 0 and 1 to represent the percentage of deuterium in solvent
      float m_DeuteriumExchangeParameter;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AnalyzeSas();

      //! @brief Clone function
      //! @return pointer to new AnalyzeAtomDistanceHeatmap
      AnalyzeSas *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return gives the string to append to the the end of a filename to identify this analysis
      const std::string &GetOutFilePostfix() const;

      //! @brief gives the experimental data file name
      //! @return the filename string
      const std::string &GetExperimentalDataFileName() const;

      //! @brief gives the sasa data file name
      //! @return the filename string
      const std::string &GetSasaDataFileName() const;

      //! @brief gives the computed data file name
      //! @return the filename string
      const std::string &GetComputedDataFileName() const;

      //! @brief gives the excluded volume parameter
      //! @return the c1 parameter
      const float &GetExcludedVolumeParameter() const;

      //! @brief gives the hydration shell thickness parameter
      //! @return the c2 parameter
      const float &GetHydrationShellParameter() const;

      //! @brief gives maximum dimension from the p(r) curve of the experimental data
      //! @return dmax
      const float &GetMaximumDimension() const;

      const float &GetC1Min() const;
      const float &GetC1Max() const;
      const float &GetC2Min() const;
      const float &GetC2Max() const;
      const float &GetC1StepSize() const;
      const float &GetC2StepSize() const;
      const float &GetDeuteriumExchangeParameter() const;

      //! @brief gives name of the score to be used to compare two saxs profiles
      //! @return Score type
      const std::string &GetScoreTypeName() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string operator()( const assemble::ProteinEnsemble &ENSEMBLE) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief hook that derived classes can override, called after TryRead successfully reads a serializer containing only initializer info (no data variables)
      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM);

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

    }; // class AnalyzeSas

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ANALYZE_SAS_H_
