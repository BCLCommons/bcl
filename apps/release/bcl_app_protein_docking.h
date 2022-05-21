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

#ifndef BCL_APP_PROTEIN_DOCKING_H_
#define BCL_APP_PROTEIN_DOCKING_H_

// include the namespace header
#include "app/bcl_app.h"

// include other forward headers - sorted alphabetically
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @class ProteinDocking
    //! @brief This is an application for predicting structure of IMP oligomers given the individual tertiary structures
    //! of their subunits. The input structures to this algorithm are a structure of the “receptor” and a structure of
    //! the “ligand”, both oriented in the membrane where the z-axis is aligned with the membrane normal using the PPM
    //! server (Lomize et al., 2012) separately. Generation of docking candidates begins with a random rotation of the
    //! ligand about the z-axis and translation of the ligand in the membrane to create a glancing contact with the receptor.
    //! The ligand is then randomly moved with respect to the receptor using a Monte Carlo search. Translation along the
    //! z-axis is limited to no more than 5.4 Å and the step size of tilt angle from the z-axis is limited to no more than
    //! 0.05 radians (or ~2.8 degree). The baseline scoring function of this algorithm consists of a clashing term that
    //! represent van der Waals repulsion, a residue pair contact potential term for interface interaction, and a radius
    //! of gyration term that favors compact packing between the two docking partners. This algorithm was designed to be
    //! able to use various restraints derived from experimental and computational approaches. In the current study,
    //! we tested the idea of using predicted interface residues and their predicted WCNs as restraints for scoring docking
    //! solutions on a set of 16 α-helical IMP oligomers.
    //!
    //! @author lib14
    //! @motified April 20th, 2018
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinDocking :
      public Interface
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const ApplicationType ProteinDocking_Instance;

    private:

      //! flag for the receptor PDB file
      util::ShPtr< command::FlagStatic> m_FlagReceptor;

      //! flag for the ligand PDB file
      util::ShPtr< command::FlagStatic> m_FlagLigand;

      //! flag for the native structure of the complex
      util::ShPtr< command::FlagStatic> m_FlagNative;

      //! flag for the membrane for the complex
      util::ShPtr< command::FlagStatic> m_FlagMembrane;

      //! flag for the file that specifies the optimizer
      util::ShPtr< command::FlagStatic> m_FlagOptimizer;

      //! flag for number of models
      util::ShPtr< command::FlagStatic> m_FlagNumberModels;

      //! flag for packing density prediction filename
      util::ShPtr< command::FlagStatic> m_FlagPackingDensity;

      //! flag for the output file prefix
      util::ShPtr< command::FlagStatic> m_FlagOutputPrefix;

      //! flag for specifying symmetric docking

    public:

      //! @brief default constructor
      ProteinDocking();

      //! @brief returns a pointer to a copy of this object
      ProteinDocking *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief a detailed description of this application
      std::string GetDescription() const;

      //! @brief initializes the command object for that executable
      //! @return initialized command object
      util::ShPtr< command::Command> InitializeCommand() const;

    protected:

    ////////////////
    // operations //
    ////////////////

      //! @brief the main function of this application
      //! @return exit code - 0 for success
      int Main() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief reads in the native protein structure and add it to ProteinModelData of MODEL
      //! @param MODEL protein model to be pre-processed
      //! @return MODEL pre-process protein model
      void PreProcess( assemble::ProteinModel &MODEL) const;

      //! @brief
      //! @param
      //!
      void PreDock( const assemble::ProteinModel &RECEPTOR, assemble::ProteinModel &LIGAND) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief reads members from a given input stream
      //! @param ISTREAM input stream to read members from
      //! @return input stream which members were read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief writes members into a given output stream
      //! @param OSTREAM output stream to write members into
      //! @param INDENT number of indentations
      //! @return output stream into which members were written
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // end of class ProteinDocking
  } // namespace app
} // namespace bcl

#endif // BCL_APP_PROTEIN_DOCKING_H_
