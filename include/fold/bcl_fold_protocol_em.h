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

#ifndef BCL_FOLD_PROTOCOL_EM_H_
#define BCL_FOLD_PROTOCOL_EM_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_protocol_interface.h"
#include "density/bcl_density_map.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_enum.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProtocolEM
    //! @brief protocol for using em density maps
    //! @details TODO: add an detailed description to this class
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Nov 18, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProtocolEM :
      public ProtocolInterface
    {

    public:

    //////////
    // data //
    //////////

      mutable Score e_BodyRestraintScore;                    //!< body restraint score
      mutable Score e_BodyConnectivityDensityScore;          //!< body connectivity score
      mutable storage::Set< Score> m_DensityAgreementScores; //!< density agreement scores

      //! density map associate with the density agreement scores
      /* util::ShPtr< density::Map> m_AgreementScoreDensityMap; */
      density::Map m_AgreementScoreDensityMap;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProtocolEM();

    public:

      //! @brief Clone function
      //! @return pointer to new ProtocolEM
      ProtocolEM *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief access to only instance
      //! @return reference to only instance
      static ProtocolEM &GetInstance();

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ///////////
    // flags //
    ///////////

      //! @brief returns all flags that are specialized for this protocol
      //! @return all flags that are specialized for this protocol
      const util::ShPtrVector< command::FlagInterface> &GetAllFlags() const;

      //! @brief return command line flag for using body restraints
      //! @return command line flag for using body restraints
      static util::ShPtr< command::FlagInterface> &GetFlagBodyRestraint();

      //! @brief return command line flag for writing body assignments when using EMFold
      //! @return command line flag for writing body assignments when using EMFold
      static util::ShPtr< command::FlagInterface> &GetFlagPrintBodyAssignment();

      //! @brief return command line flag for using density agreement score
      //! @return command line flag for using density agreement score
      static util::ShPtr< command::FlagInterface> &GetFlagScoreDensityAgreement();

      //! @brief return command line flag for using density connectivity score
      //! @return command line flag for using density connectivity score
      static util::ShPtr< command::FlagInterface> &GetFlagScoreDensityConnectivity();

      //! @brief return command line flag for using em fold refinement protocol
      //! @return command line flag for using em fold refinement protocol
      static util::ShPtr< command::FlagInterface> &GetFlagEMRefinement();

    ////////////////
    // operations //
    ////////////////

      //! @brief modifies the start model
      //! @param START_MODEL Start model to be modified
      void ModifyStartModel( assemble::ProteinModel &START_MODEL) const;

      //! @brief initialize the scores and add them to Scores enumerator
      void InitializeScores();

      //! @brief modify the score weight set
      //! @param SCORE_WEIGHT_SET Score weight set
      void ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const;

      //! @brief modify the terminate object
      //! @param CRITERION which will be modified by protocols
      //! @param STAGE Stage in which this terminate will be used
      void ModifyCriterion
      (
        opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
        const mc::Stage &STAGE
      ) const;

      //! @brief modify the printer object
      //! @param PRINTER which will be modified by protocols
      //! @param STAGE Stage in which this terminate will be used
      void ModifyPrinter
      (
        mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
        const mc::Stage &STAGE
      ) const;

      //! @brief modify the pdb factory object
      //! @param FACTORY pdb factory to be modified
      //! @param STAGE Stage in which this terminate will be used
      void ModifyFactory
      (
        util::ShPtr< pdb::Factory> &FACTORY,
        const mc::Stage &STAGE
      ) const;

      //! @brief initialize the mutates and add them to Mutates enumerator
      void InitializeMutates();

      //! @brief modify the mutate tree used
      //! @param MUTATE_TREE MutateTree to be modified
      void ModifyMutateTree( MutateTree &MUTATE_TREE) const;

      //! @brief get the mutate tree associated with this protocol
      //! @return the mutate tree associated with this protocol
      util::ShPtr< MutateTree> GetMutateTree() const;

      //! @brief merges this protocol's mutate tree into given mutate tree
      //! @param MUTATE_TREE tree into which to merge this protocol's tree
      void MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////
    // readme //
    ////////////

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

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

    public:

      //! @brief creates the mutates for adding SSEs into the Protein model into body restraints
      static void InitializeBodyAddMutates();

      //! @brief return body restraint
      //! @return body restraint
      static util::ShPtr< restraint::Body> GetBodyRestraint();

      //! @brief GetScoresDensityAgreement creates the score objects for density agreement scores
      util::ShPtrVector< score::ProteinModel> GetScoresDensityAgreement();

      //! @brief GetScoreBodyConnectivityDensity creates the score object for density connectivity score
      static util::ShPtr< score::ProteinModel> GetScoreBodyConnectivityDensity();

      //! @brief GetBodyRestraintScore creates the scoring object for use with the body restraint
      //! @return creates the scoring object for use with the body restraint
      static util::ShPtr< score::ProteinModel> GetScoreBodyRestraint();

      //! @brief singleton function for creating static body restraints to be used during folding
      static
      const util::ShPtr< util::ShPtrVector< restraint::Body> > &GetBodyRestraints();

    }; // class ProtocolEM

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_PROTOCOL_EM_H_
