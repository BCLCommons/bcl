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

#ifndef BCL_RESTRAINT_ANALYZE_ACCESSIBILITY_CHANGE_H_
#define BCL_RESTRAINT_ANALYZE_ACCESSIBILITY_CHANGE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_accessibility_profile.h"
#include "assemble/bcl_assemble_aa_exposure_interface.h"
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeAccessibilityChange
    //! @brief Analyzes accessibility changes between two different ensembles of models compared to experimental data
    //! @details Takes two ensembles of models. Calculates the mean, stddev, zscore of exposure change between them
    //!          for residues that have experimental data. The experimental data is provided by an input file. Cutoffs
    //!          (min and max) can be given for a residue to be outputted as "of interest". Outputs a text file with
    //!          these numbers and the information about the residues involved. Also, outputs a pymol script file
    //!          which colors the backbone of residues according to the experimental accessibility change and colors
    //!          the side chains of residues according to the calculated accessibility change.
    //!
    //! @see @link example_restraint_analyze_accessibility_change.cpp @endlink
    //! @author alexanns
    //! @date Nov 8, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzeAccessibilityChange :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the postfix that will be appended to the filename that will hold this analysis
      std::string m_OutFilePostFix;

      //! the ensemble of models in the starting state
      assemble::ProteinEnsemble m_StartEnsemble;

      //! the method for calculating exposure
      util::ShPtr< assemble::AAExposureInterface> m_ExposureMethod;

      //! the exposure profiles of interest
      storage::List< AccessibilityProfile> m_ExperimentalExposures;

      //! minimum mean exposure change necessary for residue to be considered
      double m_MeanMinCutoff;

      //! maximum mean exposure change allowed for residue to be considered
      double m_MeanMaxCutoff;

      //! minimum mean zscore change necessary for residue to be considered
      double m_ZScoreMinCutoff;

      //! maximum mean zscore change allowed for residue to be considered
      double m_ZScoreMaxCutoff;

      //! filename for the pymol script that will be output
      std::string m_PymolOuputFilename;

      //! index of the model that should be used as the representative for visualization in pymol
      size_t m_EnsembleRepresentativeIndex;

      //! true if the representative (and index) should be from the start ensemble - otherwise indicates from end
      bool m_EnsembleRepresentativeFromStartEnsemble;

      //! minimum value of the color gradient used for visualization in pymol
      double m_GradientMin;

      //! maximum value of the color gradient used for visualization in pymol
      double m_GradientMax;

      //! true if the experimental data and exposure calculation are directly related - otherwise inversely related
      bool m_DirectRelation;

      //! true if the user wants to explicitly set the calculated exposure ranges - otherwise automatically determined
      bool m_SetNCRange;

      //! the minimum range of calculated exposure as set by the user
      double m_NCRangeMin;

      //! the maximum range of calculated exposure as set by the user
      double m_NCRangeMax;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AnalyzeAccessibilityChange();

      //! @brief Clone function
      //! @return pointer to new AnalyzeAccessibilityChange
      AnalyzeAccessibilityChange *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return gives the string to append to the the end of a filename to identify this analysis
      const std::string &GetOutFilePostfix() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string operator()( const assemble::ProteinEnsemble &ENSEMBLE) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    public:

      //! @brief returns dummy name for ensemble
      //! @param ENSEMBLE for which a name will be created
      //! @return string which is the dummy name of ensemble
      static std::string EnsembleAsFilename
      (
        const assemble::ProteinEnsemble &ENSEMBLE
      );

      //! @brief create ensemble from filename
      //! @param ENSEMBLE ensemble to set
      //! @param NAME string name of file which will be used to create the ensemble
      //! @param ERR_STREAM stream to write out errors to
      //! @return ensemble created from the filename
      static bool
      EnsembleFromFilename( assemble::ProteinEnsemble &ENSEMBLE, const std::string &NAME, std::ostream &ERR_STREAM);

      //! @brief returns dummy name for exposure data file
      //! @return string which could be the name of the file the data comes from
      static std::string ExposureDataAsFilename
      (
        const storage::List< AccessibilityProfile> &DATA
      );

      //! @brief reads in exposure data from a file given the filename
      //! @param PROFILES a list of accessibility profiles to set
      //! @param NAME the name of the file the data will be read from
      //! @param ERR_STREAM the stream any error will be written to
      //! @return true on success
      static bool ExposureDataFromFilename
      (
        storage::List< AccessibilityProfile> &PROFILES,
        const std::string &NAME,
        std::ostream &ERR_STREAM
      );

    protected:

      //! @brief determines the range of the calculated exposures and experimental exposures
      //! @param MEAN_DIFF_STATS to hold, over all profiles, the average of the average exposure difference calculated
      //!                        across all the residues in the profiles
      //! @return vector nd with two ranges, one for calculate the other for experimental exposures
      storage::VectorND< 2, math::Range< double> > GetCalculatedAndExperimentalRanges
      (
        const assemble::ProteinEnsemble &ENSEMBLE,
        math::RunningAverageSD< double> &MEAN_DIFF_STATS
      ) const;

    }; // class AnalyzeAccessibilityChange

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ANALYZE_ACCESSIBILITY_CHANGE_H_
