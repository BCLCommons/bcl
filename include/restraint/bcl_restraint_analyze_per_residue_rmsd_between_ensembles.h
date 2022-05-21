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

#ifndef BCL_RESTRAINT_ANALYZE_PER_RESIDUE_RMSD_BETWEEN_ENSEMBLES_H_
#define BCL_RESTRAINT_ANALYZE_PER_RESIDUE_RMSD_BETWEEN_ENSEMBLES_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "quality/bcl_quality.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzePerResidueRMSDBetweenEnsembles
    //! @brief Calculates the per residue RMSD between two sets of ensembles. CA atoms are used in the calculation
    //! @details Each
    //! model in each ensemble is pairwise superimposed according to the provided superimposition method, then the
    //! per residue rmsd is calculated along the sequence. These rmsds are averaged over all pairwise RMSDs between
    //! the ensembles. Also outputs a python file that can be opened with pymol in order to vizualize the regions
    //! of the provided template structure that have large per residue RMSDs. This is depicted as thickness using
    //! the putty representation in pymol.
    //!
    //! @see @link example_restraint_analyze_per_residue_rmsd_between_ensembles.cpp @endlink
    //! @author alexanns
    //! @date Feb 22, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzePerResidueRMSDBetweenEnsembles :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the postfix that will be appended to the filename that will hold this analysis
      std::string m_OutFilePostFix;

      //! the method to use for superimposing the structures onto the template
      quality::SuperimposeMeasure m_SuperimposeMeasure;

      //! the model that the other models will be superimposed onto
      assemble::ProteinModel m_TemplateModel;

      quality::Measure m_QualityMeasure;
  
      //! filename for the pymol script that will be output
      std::string m_PymolOuputFilename;

      //! the ensemble of models in the starting state
      assemble::ProteinEnsemble m_StartEnsemble;

      //! subtract the variability from the provided ensemble from the variability between ensembles, if true
      bool m_NormalizeByInternalRMSD;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AnalyzePerResidueRMSDBetweenEnsembles();

      //! @brief Clone function
      //! @return pointer to new AnalyzePerResidueRMSDBetweenEnsembles
      AnalyzePerResidueRMSDBetweenEnsembles *Clone() const;

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

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

    }; // class AnalyzePerResidueRMSDBetweenEnsembles

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ANALYZE_PER_RESIDUE_RMSD_BETWEEN_ENSEMBLES_H_ 
