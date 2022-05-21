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

#ifndef BCL_RESTRAINT_ANALYZE_PER_RESIDUE_RMSD_H_
#define BCL_RESTRAINT_ANALYZE_PER_RESIDUE_RMSD_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "quality/bcl_quality_measures.h"
#include "quality/bcl_quality_superimpose_measures.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzePerResidueRMSD
    //! @brief Calculates the CA RMSD per residue pairwise for each model in a provided ensemble.
    //! @details For an ensemble of models, superimposes them onto a template structure using the provided superimpose
    //! method.
    //! The rmsd of each residue between the ensemble of models is then calculated. Also outputs a python script that
    //!  can be used with pymol to visualize the rmsd over the structure by using the putty representation. Residues
    //!  with larger rmsd will be thicker. The template model is used in the pymol script.
    //!
    //! @see @link example_restraint_analyze_per_residue_rmsd.cpp @endlink
    //! @author alexanns
    //! @date Feb 20, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzePerResidueRMSD :
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

      //! quality measure to use
      quality::Measure m_QualityMeasure;

      //! filename for the pymol script that will be output
      std::string m_PymolOuputFilename;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AnalyzePerResidueRMSD();

      //! @brief Clone function
      //! @return pointer to new AnalyzePerResidueRMSD
      AnalyzePerResidueRMSD *Clone() const;

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

    ///////////////
    // operators //
    ///////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string operator()( const assemble::ProteinEnsemble &ENSEMBLE) const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    public:

      //! @brief converts a protein model into a string by returning the pdb filename stored in the protein model data
      //! @param MODEL the model that will be converted into a string
      //! @return the string of the pdb filename stored in the protein model data
      static std::string ProteinModelAsString( const assemble::ProteinModel &MODEL);

      //! @brief converts a string into a protein model by assuming the string is a pdb filename
      //! @param MODEL protein model to setup
      //! @param NAME string which is the pdb filename
      //! @param ERR_STREAM stream to write out erros to
      //! @return protein model created from the pdb file
      static bool ProteinModelFromString( assemble::ProteinModel &MODEL, const std::string &NAME, std::ostream &ERR_STREAM);

      //! @brief gives the CA atom coordinates for each residue for each model in an ensemble
      //!        Each model in the ensemble is superimposed onto a templete model according to the provided
      //!        superimposition method.
      //! @param ENSEMBLE the models whose ca coordinate will be gotten for each residue
      //! @param SUPERIMPOSE_MEASURE the method of superimposition used for the models of the ensemble onto the template
      //! @param TEMPLATE_MODEL the model that the ensemble models will be superimposed on
      //! @return map with a locator for each atom and the vector of associated coordinates coming from the models in
      //!         the ensemble
      static storage::Map< assemble::LocatorAtom, storage::Vector< linal::Vector3D> > GetAtomCoordinates
      (
        const assemble::ProteinEnsemble &ENSEMBLE,
        const quality::SuperimposeMeasure &SUPERIMPOSE_MEASURE,
        const assemble::ProteinModel &TEMPLATE_MODEL
      );

    }; // class AnalyzePerResidueRMSD

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ANALYZE_PER_RESIDUE_RMSD_H_
