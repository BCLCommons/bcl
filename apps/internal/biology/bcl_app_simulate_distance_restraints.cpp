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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_collector_aa_specified.h"
#include "assemble/bcl_assemble_locator_atom.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "math/bcl_math_histogram.h"
#include "nmr/bcl_nmr_star_noe_handler.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "random/bcl_random_histogram_1d_distribution.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "restraint/bcl_restraint_handler_atom_distance_assigned.h"
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes

namespace bcl
{
  namespace app
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintSimulateDistances
    //! @brief app is for creating an AtomDistanceAssigned restraint file.
    //!
    //! @author alexanns, teixeipl
    //! @date 10/18/2008
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintSimulateDistances :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! flag for pdb file from which the restraints will be simulated
      util::ShPtr< command::FlagInterface> m_PDBFilename;

      //! flag for specifying that distance restraints should be simulated
      util::ShPtr< command::FlagInterface> m_SimulateDistanceRestraint;

      //! flag for specifying that simulated distance restraints will not consider undefined AAs in the given pdb file
      util::ShPtr< command::FlagInterface> m_SkipUndefinedAAs;

      //! flag for specifying that NMR distance restraints should be simulated
      util::ShPtr< command::FlagStatic> m_SimulateNMRDistanceRestraint;
      util::ShPtr< command::ParameterInterface> m_AaDistanceParam;
      util::ShPtr< command::ParameterInterface> m_NMRILVOnlyRestraints;

      //! flag for the specifiying an exact number of desired restraints
      util::ShPtr< command::FlagInterface> m_NumberDesiredRestraints;

      //! flag for specifying the upperlimit a distance restraint can have
      util::ShPtr< command::FlagStatic> m_Limit;
      util::ShPtr< command::ParameterInterface> m_UpperLimit;
      util::ShPtr< command::ParameterInterface> m_LowerLimit;

      //! flag specifying whether restraints should be ignored if within a given minimum sequence separation, default 0
      util::ShPtr< command::FlagStatic> m_MinimumSeparation;

      //! flag for specifying a maximum neighbor count limit for an atom involved in a restraint
      util::ShPtr< command::FlagInterface> m_NeighborCountLimit;

      //! flag for specifying the output file path and name
      util::ShPtr< command::FlagInterface> m_OutputFile;

      //! flag that only residues between specified lists should be considered for restraints
      util::ShPtr< command::FlagStatic> m_ConsiderExplicitResidues;
      util::ShPtr< command::ParameterInterface> m_ExplicitResidueListFilenameA;
      util::ShPtr< command::ParameterInterface> m_ExplicitResidueListFilenameB;

      //! flag specifying that restraints involving residues changing the most between two states should be gotten
      util::ShPtr< command::FlagStatic> m_OrderByDistanceChange;
      util::ShPtr< command::ParameterInterface> m_EndStatePDBFilename; //< the pdb with the second state

      //! flag specifying restraints should be written in rosetta 2.0 format
      util::ShPtr< command::FlagStatic> m_WriteRosettaFormattedRestraints;

      //! flag specifying that a random amount should be added to the first side chain atom distance between two
      //! residues. This can be used to simulate the uncertainty in a epr distance measurement.
      util::ShPtr< command::FlagStatic> m_AddDistributionBiasToDistance;
      //! This is the filename containing the histogram distribution which will be used to come up with the random
      //! values
      util::ShPtr< command::ParameterInterface> m_DistributionBiasFilename;

      //! flag specifying restraints should be written in rosetta 3.0 format
      util::ShPtr< command::FlagStatic> m_WriteRosettaMiniFormattedRestraints;

      //! flag specifying restraints should be written in modified CASP RR format
      util::ShPtr< command::FlagStatic> m_WriteModifiedCASPFormattedRestraints;

      //! flag specifying that a list of desired restraints will be provided
      util::ShPtr< command::FlagStatic> m_RestraintListProvided;
      util::ShPtr< command::ParameterInterface> m_RestraintListFilename;     //< filename containing restraints
      util::ShPtr< command::ParameterInterface> m_RestraintListChainAColumn; //<
      util::ShPtr< command::ParameterInterface> m_RestraintListResiAColumn;  //<
      util::ShPtr< command::ParameterInterface> m_RestraintListChainBColumn; //<
      util::ShPtr< command::ParameterInterface> m_RestraintListResiBColumn;  //<

      //! flag specifying that the number of restraints should be a fraction of the residues in the protein
      mutable util::ShPtr< command::FlagStatic> m_FractionDesiredRestraints;
      util::ShPtr< command::ParameterInterface> m_FractionDesiredRestraintsNumber;

      //! map to store the residue identities for all input residues
      mutable storage::Map< storage::Pair< char, int>, char> m_ResidueMap;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RestraintSimulateDistances();

    public:

      //! @brief Clone function
      //! @return pointer to new Quality
      RestraintSimulateDistances *Clone() const
      {
        return new RestraintSimulateDistances( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const
      {
        return storage::Vector< std::string>( size_t( 1), "SimulateDistanceRestraints");
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief GenerateProteinModel creates and returns a protein model based on a command flag
      //! @param PDB_STRING is the string which gives the pdb to create the protein model from
      //! @return ProteinModel which was created from the pdb given by "PDB_STRING"
      assemble::ProteinModel GenerateProteinModel( const std::string &PDB_STRING) const;

      //! @brief WriteDistanceRestraintsRosettaFormat writes a list of restraint information in rosetta format
      //! @param RESTRAINTS the list of restraint information
      //! @param OSTREAM the stream to which the restraint will be written
      //! @return std::ostream which was passed as parameter
      std::ostream &WriteDistanceRestraintsRosettaFormat
      (
        const storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        > &RESTRAINTS,
        std::ostream &OSTREAM
      ) const;

      //! @brief GetAccessibleResidues determines the residues from a protein model that meet a neighbor count cutoff
      //! @param MODEL the protein model from which the residues will be taken
      //! @return vector of amino acids which meet a neighbor count cutoff
      const util::SiPtrVector< const biol::AABase> GetAccessibleResidues( const assemble::ProteinModel &MODEL) const;

      //! @brief GetAllRestraintsInfo creates restraint information all pairs of a list of residues
      //! @param RESIDUES the residues for which pairwise restraint information will be created
      //! @return list which has the information for all pairwise restraints in RESIDUES
      const storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > GetAllRestraintsInfo( const util::SiPtrVector< const biol::AABase> &RESIDUES) const;

      //! @brief RandomizeRestraintOrder takes a vector of restraint information and randomizes its order
      //! @param RESTRAINTS the vector of restraint information whose order will be randomized
      //! @return vector with the restraints in a random order
      storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > RandomizeRestraintOrder
      (
        const storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        > &RESTRAINTS
      ) const;

      //! @brief OrderRestraintByDistanceChange calculates the distance change from start to end state and orders the
      //!        restraints so that those with the largest distance change are first
      //! @param RESTRAINTS the restraints that will be ordered
      //! @param END_STATE the other state of the protein for which distances will be measured and will be used to get
      //!        distance changes
      //! @return vector of restraints ordered so that the first restraint has the largest distance change in going from
      //!         the starting state to the ending state
      storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > OrderRestraintByDistanceChange
      (
        const storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        > &RESTRAINTS,
        const assemble::ProteinModel &START_STATE,
        const assemble::ProteinModel &END_STATE
      ) const;

      //! @brief WriteNumberOfDesiredRestraints outputs the desired number of restraints in the desired format to the
      //!         desired filename
      //! @param RESTRAINTS the restraints some of which will be output
      //! @param NUMBER_TO_OUTPUT the number of restraints that will be output
      void WriteNumberOfDesiredRestraints
      (
        const storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        > &RESTRAINTS,
        const size_t NUMBER_TO_OUTPUT
      ) const;

      //! @brief GetExplicitlyAllowedRestraints gives the restraint information for restraints based on the residue lists
      //!        passed over the command line. The restraints have to involve residues that are also in the provide list.
      //!        So restraint information is created for all residue pairs between the two command line lists, but then
      //!        restraints are removed that involve residues not contained in the argument parameter residue list
      //! @param RESIDUES the list of residues that can be involved in restraints that the restraints are filtered
      //!        according to
      //! @return vector of restraint information
      storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > GetExplicitlyAllowedRestraints
      (
        const util::SiPtrVector< const biol::AABase> &RESIDUES
      ) const;

      //! @brief GetExactRestraints gives the restraint information for restraints based on the list of restraints
      //!        passed over the command line. The list contains all pairwise restraints that will be created.
      //!        Restraints are removed that involve residues not contained in the argument parameter residue list
      //!        (which is usually the list of exposed residues).
      //!        The format of the file should be one restraint per line as :
      //!        <chain_col> <resi_seq_id_col> <chain_col> <resi_seq_id_col>
      //!        an example would be :
      //!        'A' 23 'A' 46
      //!        'B' 19 'C' 10
      //! @param RESIDUES the list of residues that can be involved in restraints that the restraints are filtered
      //!        according to
      //! @return vector of restraint information
      storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > GetExactRestraints
      (
        const util::SiPtrVector< const biol::AABase> &RESIDUES
      ) const;

      //! @brief GetDistanceUpperLowerBounds gives the distance between the first two coordinates as well as the upper
      //!        and lower bound
      //! @param COORDINATES_A the first coordinate involved in the restraint
      //! @param COORDINATES_B the second coordinate involved in the restraint
      //! @return VectorND< 3> with the distance, upper bound, and lower bound, respectively, of the restraint
      storage::VectorND< 3, double> GetDistanceUpperLowerBounds
      (
        const linal::Vector3D &COORDINATES_A, const linal::Vector3D &COORDINATES_B
      ) const;

      //! @brief GetRestraintInformation gets the information necessary for a restraint from two residues and the
      //!        distance, upper, and lower bounds
      //! @param RESI_A the first residue involved in the restraint
      //! @param RESI_B the second residue involved in the restraint
      //! @param DISTANCE_UPPER_LOWER_BOUND the distance, upper bound, lower bound, respectively of the restraint
      //! @return the information about the restraint
      storage::Pair
      <
        storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
      > GetRestraintInformation
      (
        const biol::AABase &RESI_A, const biol::AABase &RESI_B,
        const storage::VectorND< 3, double> &DISTANCE_UPPER_LOWER_BOUND
      ) const;

      //! @brief GetDistributionBiasAmount provides a value based on the distribution provided over the command line
      //! @return double which is the random value obtained from the distribution
      double GetDistributionBiasAmount() const;

      //! @brief WriteDistanceRestraintsModifiedCASPFormat writes a list of restraint information in modified CASP format
      //! @param RESTRAINTS the list of restraint information
      //! @param OSTREAM the stream to which the restraint will be written
      //! @return std::ostream which was passed as parameter
      std::ostream &WriteDistanceRestraintsModifiedCASPFormat
      (
        const storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        > &RESTRAINTS,
        std::ostream &OSTREAM
      ) const;

      //! @brief WriteDistanceRestraintsRosettaMiniFormat writes a list of restraint information in rosetta mini format
      //! @param RESTRAINTS the list of restraint information
      //! @param OSTREAM the stream to which the restraint will be written
      //! @return std::ostream which was passed as parameter
      std::ostream &WriteDistanceRestraintsRosettaMiniFormat
      (
        const storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        > &RESTRAINTS,
        std::ostream &OSTREAM
      ) const;

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
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief creates a shptr vector of Atom Distance Restraints to be written
      //! @param PROTEIN_MODEL where list of residues and restraint distances come from, must be protonated
      //!        and AAComplete
      //! @param NUMBER_DESIRED_RESTRAINTS how many restraints should be generated
      //! @param UPPER_LIMIT_DEVIATION the amount which should be added to the distance found to find the upper bound
      //! @param AA_DISTANCE the smallest sequence distance two aa's can be in sequence for a restraint to be made
      //! @param DISTANCE_RANGE all distances must be in this range to be considered
      //! @return a shptr vector of Atom Distances
      util::ShPtrVector< restraint::AtomDistance> CreateNOERestraints
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        const double NUMBER_DESIRED_RESTRAINTS,
        const double UPPER_LIMIT_DEVIATION,
        const size_t AA_DISTANCE,
        const math::Range< double> &DISTANCE_RANGE
      ) const;

      //! @brief takes a protein model and returns all the atom pairs that have distances within the range
      //! @param PROTEIN_MODEL where list of residues and restraint distances come from, must be protonated
      //!        and AAComplete
      //! @param AA_DISTANCE the smallest sequence distance two aa's can be in sequence for a restraint to be made
      //! @param DISTANCE_RANGE all distances must be in this range to be considered
      //! @return all the atom pairs that have distances within the range
      storage::List
      <
        storage::VectorND< 2, storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> > >
      > GetAllProtonPairs
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        const size_t AA_DISTANCE,
        const math::Range< double> &DISTANCE_RANGE
      ) const;

      //! @brief returns whether this atom type is valid for ILV-labelling (methyl in ILV or backbone H otherwise)
      //! @param AA_TYPE AA type this atom is associatied with
      //! @param ATOM_TYPE atom type of the atom
      //! @return whether this atom type is valid for ILV-labelling (methyl in ILV or backbone H otherwise)
      static bool IsValidILVAtom( const biol::AAType &AA_TYPE, const biol::AtomType &ATOM_TYPE);

    private:

      static const ApplicationType RestraintSimulateDistances_Instance;

    }; // RestraintSimulateDistances

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> RestraintSimulateDistances::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add flag
      sp_cmd->AddFlag( m_PDBFilename);

      // add flag
      sp_cmd->AddFlag( m_SimulateDistanceRestraint);

      // add flag
      sp_cmd->AddFlag( m_SkipUndefinedAAs);

      // add flag
      sp_cmd->AddFlag( m_SimulateNMRDistanceRestraint);

      // add flag
      sp_cmd->AddFlag( m_NumberDesiredRestraints);

      // add flag
      m_FractionDesiredRestraints->PushBack( m_FractionDesiredRestraintsNumber);
      sp_cmd->AddFlag( m_FractionDesiredRestraints);

      // add flag
      sp_cmd->AddFlag( m_Limit);

      // add flag
      sp_cmd->AddFlag( m_MinimumSeparation);

      // add flag
      sp_cmd->AddFlag( m_NeighborCountLimit);

      // add flag
      sp_cmd->AddFlag( m_OutputFile);

      // add flag for explicit residues
      sp_cmd->AddFlag( m_ConsiderExplicitResidues);

      // add flag flag specifying restraints involving residues changing the most between two states should be gotten
      sp_cmd->AddFlag( m_OrderByDistanceChange);

      // add flag specifying restraints should be written in rosetta 2.0 format
      sp_cmd->AddFlag( m_WriteRosettaFormattedRestraints);

      // flag specifying that a random amount should be added to the first side chain atom distance between two aas
      sp_cmd->AddFlag( m_AddDistributionBiasToDistance);

      // flag specifying restraints should be written in rosetta 3.0 format
      sp_cmd->AddFlag( m_WriteRosettaMiniFormattedRestraints);

      //! flag specifying restraints should be written in modified CASP RR format
      sp_cmd->AddFlag( m_WriteModifiedCASPFormattedRestraints);

      // flag specifying that a list of desired restraints will be provided
      sp_cmd->AddFlag( m_RestraintListProvided);

      // add flag
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int RestraintSimulateDistances::Main() const
    {
      // create ProteinModel "protein_model" and initialize with protein of "pdb_filename"
      const assemble::ProteinModel
        protein_model( GenerateProteinModel( m_PDBFilename->GetFirstParameter()->GetValue()));

      // default number of restraints
      size_t number_restraints( m_NumberDesiredRestraints->GetFirstParameter()->GetNumericalValue< size_t>());

      // true if user specified fraction of residues for number of restraints
      if( m_FractionDesiredRestraints->GetFlag())
      {
        const size_t number_resi( protein_model.GetNumberAAs());
        const double fraction_desired( m_FractionDesiredRestraintsNumber->GetNumericalValue< double>());
        number_restraints = fraction_desired * double( number_resi);
        BCL_MessageStd( "fraction of restraints: " + util::Format()( fraction_desired));
        BCL_MessageStd( "number residues " + util::Format()( number_resi));
        BCL_MessageStd( "number restraints " + util::Format()( number_restraints));
      }

      // true if star file formatted restraints should be written
      if( m_SimulateNMRDistanceRestraint->GetFlag())
      {
        io::OFStream ofstream;
        io::File::MustOpenOFStream( ofstream, m_OutputFile->GetFirstParameter()->GetValue());

        // calculate the number of restraints based on the number of AAs in the protein if the number is set to default
        if( number_restraints == std::numeric_limits< size_t>::max())
        {
          number_restraints = protein_model.GetNumberAAs();
        }

        BCL_MessageDbg( "restraint number is: " + util::Format()( number_restraints));
        nmr::StarNOEHandler handler_noe;
        handler_noe.WriteRestraints
        (
          ofstream,
          CreateNOERestraints
          (
            protein_model,
            number_restraints,
            0.5,
            m_SimulateNMRDistanceRestraint->GetFirstParameter()->GetNumericalValue< size_t>(),
            math::Range< double>( m_LowerLimit->GetNumericalValue< double>(), m_UpperLimit->GetNumericalValue< double>())
          )
        );
        io::File::CloseClearFStream( ofstream);

        return 0;
      }

      // true if restraints should be simulated
      if( m_SimulateDistanceRestraint->GetFlag())
      {
        // get the residues from "protein_model" which meet the exposure cutoff
        const util::SiPtrVector< const biol::AABase> exposed_amino_acids( GetAccessibleResidues( protein_model));

//        BCL_MessageDbg( "Test before getacessibleresitudues FUNCTOIN: " + util::Format()( protein_model.GetAminoAcids().FirstElement()->GetAtomCoordinates()));
//        BCL_MessageDbg( "TESTFIRST EXPOSED AA: " + util::Format()( exposed_amino_acids.FirstElement()->GetAtomCoordinates())); // DEBUG REMOVE LINE

        // create vector which will hold all of the restraint information
        storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        > all_restraints_info;

        // true if only certain residues should be considered for restraints
        if( m_ConsiderExplicitResidues->GetFlag())
        {
          BCL_MessageDbg( "First"); // DEBUG REMOVE LINE
          all_restraints_info = GetExplicitlyAllowedRestraints( exposed_amino_acids);
        }
        else if( m_RestraintListProvided->GetFlag())
        {
          BCL_MessageDbg( "second"); // DEBUG REMOVE LINE
          all_restraints_info = GetExactRestraints( exposed_amino_acids);
        }
        else //< get all restraints
        {
          BCL_MessageDbg( "third"); // DEBUG REMOVE LINE
          all_restraints_info = GetAllRestraintsInfo( exposed_amino_acids);
        }

        // filter out all restraints that are within minimum separation if given
        const size_t min_separation( m_MinimumSeparation->GetFirstParameter()->GetNumericalValue< size_t>());
        if( min_separation > 0)
        {
          // create vector which will hold the temporary filtered restraint information
          storage::Vector
          <
            storage::Pair
            <
              storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
            >
          > filtered_restraints_info;

          // set up iterators
          storage::Vector
          <
            storage::Pair
            <
              storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
            >
          >::const_iterator itr( all_restraints_info.Begin()), itr_end( all_restraints_info.End());

          for( ; itr < itr_end; ++itr)
          {
            // ensure that the restraint sequence separation is greater than the minimum allowed separation
            if
            (
              // Compares first amino acid position with the second amino acid position and chain
              itr->First().First().First() == itr->First().Second().First()
              && abs( itr->First().Second().Second() - itr->First().First().Second()) <= int( min_separation)
            )
            {
              continue;
            }

            filtered_restraints_info.PushBack( *itr);
          }

          // reset the restraints list using the newly filtered list
          all_restraints_info = filtered_restraints_info;
        }
        // make sure some restraints will be written
        BCL_Assert
        (
          !all_restraints_info.IsEmpty(),
          "there are no restraints that match the NC cutoff and are explicitly allowed (if the flag was used)"
        );

        // true if the restraints should be ordered by distance change
        if( m_OrderByDistanceChange->GetFlag())
        {
          all_restraints_info = OrderRestraintByDistanceChange
          (
            all_restraints_info,
            protein_model, //< start state
            GenerateProteinModel( m_EndStatePDBFilename->GetValue()) //< end state
          );
        }
        else //< randomly order the restraints
        {
          all_restraints_info = RandomizeRestraintOrder( all_restraints_info);
        }

        BCL_MessageStd( "total restraints: " + util::Format()( all_restraints_info.GetSize()));
        BCL_MessageStd( "number of desired restraints : " + util::Format()( number_restraints));
        // write out the restraints - as many as specified by command line
        WriteNumberOfDesiredRestraints( all_restraints_info, number_restraints);
      }

      //successful end
      return 0;
    }

    //! @brief GenerateProteinModel creates and returns a protein model based on a command flag
    //! @param PDB_STRING is the string which gives the pdb to create the protein model from
    //! @return ProteinModel which was created from the pdb given by "PDB_STRING"
    assemble::ProteinModel RestraintSimulateDistances::GenerateProteinModel( const std::string &PDB_STRING) const
    {
      // create "factory" to create protein model with amino acids of type AABackBone
      const pdb::Factory factory;

      // initialize write and read stream objects
      io::IFStream read;
      // open pdb file a
      io::File::MustOpenIFStream
      (
        read,
        PDB_STRING
      );

      // true is used to advice handler to ignore clashes in the structure and insert residues as they are suggested by
      // the numbering without regard for the bond information
      pdb::Handler pdb_a( read, true);
      io::File::CloseClearFStream( read);
      // create ProteinModel "protein_model_a" from "pdb_a"
      BCL_MessageTop( "building " + PDB_STRING + " from pdb information");

      // return protein model made from "pdb_a"
      return assemble::ProteinModel( factory.ProteinModelFromPDB( pdb_a));
    }

    //! @brief WriteDistanceRestraintsRosettaFormat writes a list of restraint information in rosetta format
    //! @param RESTRAINTS the list of restraint information
    //! @param OSTREAM the stream to which the restraint will be written
    //! @return std::ostream which was passed as parameter
    std::ostream &RestraintSimulateDistances::WriteDistanceRestraintsRosettaFormat
    (
      const storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > &RESTRAINTS,
      std::ostream &OSTREAM
    ) const
    {

      OSTREAM << "Rosetta++ Formatted" << '\n'; //< comment line
      OSTREAM << "Distance restraints" << '\n'; //< comment line
      OSTREAM << "THIRD COMMENT LINE" << '\n'; //< comment line
      OSTREAM << RESTRAINTS.GetSize() << '\n'; //< number of restraints to read

      // if this is anything other than a blank character than rosetta treats it as a comment line and ignores it
      const char comment_tag( ' ');

      for
      (
        storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        >::const_iterator restraint_itr( RESTRAINTS.Begin()), restraint_itr_end( RESTRAINTS.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        OSTREAM << comment_tag
        << "  " //< skip TWO spaces
        << util::Format().W( 4).R()( restraint_itr->First().First().Second()) //< first residue number
        << ' '
        << util::Format().W( 4).R()( restraint_itr->First().First().Third().GetType().GetName()) //< first atom type
        << ' '
        << util::Format().W( 4).R()( restraint_itr->First().Second().Second()) //< second residue number
        << ' '
        << util::Format().W( 4).R()( restraint_itr->First().Second().Third().GetType().GetName()) //< second atom type
        << ' '
        << util::Format().W( 10).R()( restraint_itr->Second().Second()) //< upper bound
        << ' '
        << util::Format().W( 10).R()( restraint_itr->Second().Third()) //< lower bound
        << ' '
        << util::Format().W( 10).R()( restraint_itr->Second().First()) //< actual distance
        << '\n';
      }

      return OSTREAM;
    }

    //! @brief WriteDistanceRestraintsModifiedCASPFormat writes a list of restraint information in modified CASP format
    //! @param RESTRAINTS the list of restraint information
    //! @param OSTREAM the stream to which the restraint will be written
    //! @return std::ostream which was passed as parameter
    std::ostream &RestraintSimulateDistances::WriteDistanceRestraintsModifiedCASPFormat
    (
      const storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > &RESTRAINTS,
      std::ostream &OSTREAM
    ) const
    {
      OSTREAM << "REMARK BCL Contacts" << '\n';

      // initialize counter
      size_t counter( 1);

      // iterate over the restraints
      for
      (
        storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        >::const_iterator itr( RESTRAINTS.Begin()), itr_end( RESTRAINTS.End());
        itr < itr_end;
        ++itr, ++counter
      )
      {
        // Triplet reference
        const storage::Triplet< char, int, biol::Atom> &triplet_a( itr->First()( 0));
        const storage::Triplet< char, int, biol::Atom> &triplet_b( itr->First()( 1));

        storage::Pair< char, int> key_a( triplet_a.First(), triplet_a.Second());

        // get first residue type char
        const char type_a( m_ResidueMap.Find( key_a)->second);

        // get second residue type char
        const char type_b( m_ResidueMap.Find( storage::Pair< char, int>( triplet_b.First(), triplet_b.Second()))->second);

        // create the output line with format CONTACT_N: 1P:2T A B 1.0
        // where P and T are aa types, A and B are chain ID and 1.0 is confidence
        OSTREAM << "CONTACT_"
          << util::Format()( counter)
          << std::string( ": ")
          << itr->First()( 0).Second()
          << util::Format()( type_a)
          << ":"
          << itr->First()( 1).Second() << util::Format()( type_b) << " "
          << itr->First()( 0).First() << " "
          << itr->First()( 1).First() << " 1.0\n";
      }

      // end
      return OSTREAM;
    }

    //! @brief WriteDistanceRestraintsRosettaMiniFormat writes a list of restraint information in rosetta mini format
    //! @param RESTRAINTS the list of restraint information
    //! @param OSTREAM the stream to which the restraint will be written
    //! @return std::ostream which was passed as parameter
    std::ostream &RestraintSimulateDistances::WriteDistanceRestraintsRosettaMiniFormat
    (
      const storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > &RESTRAINTS,
      std::ostream &OSTREAM
    ) const
    {
      // iterate through the vector of restraint information in order to write it to "OSTREAM"
      for
      (
        storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        >::const_iterator restraint_itr( RESTRAINTS.Begin()), restraint_itr_end( RESTRAINTS.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        OSTREAM << util::Format().W( 10).L()( "AtomPair")
        << util::Format().W( 4).R()( restraint_itr->First().First().Third().GetType().GetName()) //< first atom type
        << util::Format().W( 4).R()( restraint_itr->First().First().Second()) //< first residue number
        << util::Format().W( 4).R()( restraint_itr->First().Second().Third().GetType().GetName()) //< second atom type
        << util::Format().W( 4).R()( restraint_itr->First().Second().Second()) //< second residue number
        << "  SPLINE  " << "EPR_DISTANCE  "
        << util::Format().W( 10).R()( restraint_itr->Second().First()) //< actual distance
        << " 1.0 " //< restraint weight
        << " 0.5"  //< bin size of epr histogram
        << '\n';
      }

      return OSTREAM;
    }

    //! @brief GetAccessibleResidues determines the residues from a protein model that meet a neighbor count cutoff
    //! @param MODEL the protein model from which the residues will be taken
    //! @return vector of amino acids which meet a neighbor count cutoff
    const util::SiPtrVector< const biol::AABase>
    RestraintSimulateDistances::GetAccessibleResidues( const assemble::ProteinModel &MODEL) const
    {
      // create SiPtrVector "amino_acids" and initialize with the amino acids in the SSEs of "protein_model"
      const util::SiPtrVector< const biol::AABase> amino_acids( MODEL.GetAminoAcids());

      BCL_MessageDbg( "FIRST OF THE AMINO_ACIDS!: " + util::Format()( amino_acids.Begin()->operator *().GetAtomCoordinates())); // DEBUG REMOVE LINE

      BCL_MessageDbg
      (
        "has " + util::Format()( amino_acids.GetSize()) + " amino acids in its SSEs."
      );

      assemble::AANeighborCount neighbor_counter;

      // create AllAANeighborList "neighbor_list" initialize with the neighbors for each amino acid in "amino_acids"
      const assemble::AANeighborListContainer neighbor_list
      (
        amino_acids,
        neighbor_counter.GetDistanceCutoff(),
        neighbor_counter.GetMinimalSequenceSeparation(),
        true
      );

      // create double "max_neighbor_count" and initialize with the neighbor count limit passed over the command line
      const double max_neighbor_count( m_NeighborCountLimit->GetFirstParameter()->GetNumericalValue< double>());

      // print the neighbor count limit "max_neighbor_count"
      BCL_MessageDbg( "the neighbor count limit is " + util::Format()( max_neighbor_count));

      // create SiPtrVector "exposed_amino_acids" to hold the amino acids which have an acceptable neighbor count
      util::SiPtrVector< const biol::AABase> exposed_amino_acids;

      // get amino acids which are eligable in terms of neighbor count
      for
      (
        assemble::AANeighborListContainer::const_iterator
          itr( neighbor_list.Begin()), itr_end( neighbor_list.End());
        itr != itr_end;
        ++itr
      )
      {
        // create double "current_neighbor_count" and initialize with the number of neighbors around the amino acid
        // currently denoted by "itr"
        const double current_neighbor_count( neighbor_counter( itr->second));

        // write message
        BCL_MessageDbg
        (
          "the neighbor count of amino acid" + itr->second.GetCenterAminoAcid()->GetIdentification()
          + " is " + util::Format()( current_neighbor_count)
        );

        // true if "current_neighbor_count" is low enough
        if( current_neighbor_count < max_neighbor_count)
        {
          BCL_MessageDbg
          (
            "adding " + util::Format()( itr->second.GetCenterAminoAcid()->GetIdentification()) + " to exposed_amino_acids"
          );

          // add the amino acid currently denoted by "itr" to "exposed_amino_acids"
          exposed_amino_acids.PushBack( itr->second.GetCenterAminoAcid());
        }
      }

      // write message giving number of exposed amino acids
      BCL_MessageDbg
      (
        "The number of exposed amino acids is " + util::Format()( exposed_amino_acids.GetSize())
      );

      return exposed_amino_acids;
    }

    //! @brief GetAllRestraintsInfo creates restraint information all pairs of a list of residues outside of the minimum separation if set
    //! @param RESIDUES the residues for which pairwise restraint information will be created
    //! @return list which has the information for all pairwise restraints in RESIDUES
    const storage::Vector
    <
      storage::Pair
      <
        storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
      >
    > RestraintSimulateDistances::GetAllRestraintsInfo
    (
      const util::SiPtrVector< const biol::AABase> &RESIDUES
    ) const
    {
      // create list that will hold all of the pairwise restraints that meet the distance criteria
      storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > restraints;

      BCL_MessageDbg( "TESTFIRST in RESIDUES AA: " + util::Format()( RESIDUES.FirstElement()->GetAtomCoordinates())); // DEBUG REMOVE LINE

      size_t countera( 0); // DEBUG LINE REMOVE
      size_t counterb( 0); // DEBUG LINE REMOVE
      BCL_MessageStd( "COUNTER");
      BCL_MessageStd( "COUNTER aTRSESTINt: " + util::Format()( countera));

      // iterate through "RESIDUES" to build up "restraints"
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          itr( RESIDUES.Begin()), itr_end( RESIDUES.End());
        itr != itr_end;
        ++itr,
        ++countera
      )
      {
        // create SiPtr to AABase "amino_acid_a" initialize with the amino acid SiPtr indicated by "itr"
        const util::SiPtr< const biol::AABase> &amino_acid_a( *itr);

        // iterate through "RESIDUES" to build up "restraints"
        for
        (
          util::SiPtrVector< const biol::AABase>::const_iterator itr_b( itr + 1);
          itr_b < itr_end;
          ++itr_b,
          ++counterb
        )
        {
          // create SiPtr to AABase "amino_acid_b" initialize with the amino acid SiPtr indicated by "itr_b"
          const util::SiPtr< const biol::AABase> &amino_acid_b( *itr_b);

          BCL_MessageVrb( "COUNTER A at: " + util::Format()( countera));
          BCL_MessageVrb( "COUNTER B at: " + util::Format()( counterb));
//          BCL_MessageVrb( util::Format()( amino_acid_a))// DEBUG LINE REMOVE
//          BCL_MessageVrb( util::Format()( amino_acid_b))// DEBUG LINE REMOVE

          // if flag for skipping undefined AAs is set, skip undefined AA
          if
          (
            m_SkipUndefinedAAs->GetFlag()
            &&
            (
              !amino_acid_a->GetFirstSidechainAtom().GetCoordinates().IsDefined()
              || !amino_acid_b->GetFirstSidechainAtom().GetCoordinates().IsDefined()
            )
          )
          {
            BCL_MessageVrb( "CONTINUE FIRED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"); // DEBUG REMOVE LINE
            continue;
          }

          // get the distance between resi_a and resi_b first side chain atoms
          storage::VectorND< 3, double> distance_upper_lower_bound
          (
            GetDistanceUpperLowerBounds
            (
              amino_acid_a->GetFirstSidechainAtom().GetCoordinates(),
              amino_acid_b->GetFirstSidechainAtom().GetCoordinates()
            )
          );

          // DEBUG REMOVE LINE - Condition below seems to be superfluous, consider removing?#######################################
          // true if the distance is not defined. For example could be undefined if it is outside the boundaries
          // of distances that should be considered for restraints
          if( !util::IsDefined( distance_upper_lower_bound.First()))
          {
            // continue to next iteration since don't want to add this restraint if distance is undefined
            continue;
          }

          // continue since you don't want restraints whose distance is greater than the upper limit
          if( distance_upper_lower_bound.First() > m_UpperLimit->GetNumericalValue< double>())
          {
            continue;
          }

          // add the current restraint information to "restraints"
          restraints.PushBack( GetRestraintInformation( *amino_acid_a, *amino_acid_b, distance_upper_lower_bound));
        }
      }

      return restraints;
    }

    //! @brief RandomizeRestraintOrder takes a vector of restraint information and randomizes its order
    //! @param RESTRAINTS the vector of restraint information whose order will be randomized
    //! @return vector with the restraints in a random order
    storage::Vector
    <
      storage::Pair
      <
        storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
      >
    > RestraintSimulateDistances::RandomizeRestraintOrder
    (
      const storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > &RESTRAINTS
    ) const
    {
      // make a copy of "RESTRAINTS" so the the elements can be removed as they are selected so they aren't chosen twice
      storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > restraints_copy( RESTRAINTS);

      // create vector which will hold the restraints in a random order
      storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > randomized_restraints;

      // continue as long as "restraints_copy" still has at least one element
      while( !restraints_copy.IsEmpty())
      {
        // get a random index that will be used to get an element from "restraints_copy"
        size_t index( random::GetGlobalRandom().Random< size_t>( restraints_copy.GetSize() - 1));

        // add the random element denoted by "index" of "restraints_copy" into "randomized_restraints
        randomized_restraints.PushBack( restraints_copy( index));

        // remove the element located at "index" from "restraint_copy"
        restraints_copy.RemoveElements( index, 1);
      }

      // make sure "randomized_restraints" has all the restraints that "RESTRAINTS" does
      BCL_Assert
      (
        randomized_restraints.GetSize() == RESTRAINTS.GetSize(),
        "should be size " + util::Format()( RESTRAINTS.GetSize()) + " but is size " +
        util::Format()( randomized_restraints.GetSize())
      );

      // return the list of restraints in a random order
      return randomized_restraints;
    }

    //! @brief OrderRestraintByDistanceChange calculates the distance change from start to end state and orders the
    //!        restraints so that those with the largest distance change are first
    //! @param RESTRAINTS the restraints that will be ordered
    //! @param END_STATE the other state of the protein for which distances will be measured and will be used to get
    //!        distance changes
    //! @return vector of restraints ordered so that the first restraint has the largest distance change in going from
    //!         the starting state to the ending state
    storage::Vector
    <
      storage::Pair
      <
        storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
      >
    > RestraintSimulateDistances::OrderRestraintByDistanceChange
    (
      const storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > &RESTRAINTS,
      const assemble::ProteinModel &START_STATE,
      const assemble::ProteinModel &END_STATE
    ) const
    {
      // create multimap to hold all of the distance changes and the associated restraints
      std::multimap
      <
        double,
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > distance_changes_and_info;

      // iterate through the restraints and get the distance that occurs in the ending model
      // all of the distances for the starting model are already contained in the restraints
      for
      (
        storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        >::const_iterator restraint_itr( RESTRAINTS.Begin()), restraint_itr_end( RESTRAINTS.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        // create atom locator to locate first atom in restraint denoted by "restraint_itr"
        const assemble::LocatorAtom atom_locator_a
        (
          restraint_itr->First().First().First(), //< chain
          restraint_itr->First().First().Second(), //< seq id
          restraint_itr->First().First().Third().GetType() //< atom type
        );

        // get pointer to atom located with "atom_locator_a"
        const util::SiPtr< const biol::Atom> start_atom_a( atom_locator_a.LocateAtom( START_STATE));

        // make sure atom could be located in "START_STATE"
        BCL_Assert( start_atom_a.IsDefined(), "start_atom_a is not defined");

        // get pointer to atom located with "atom_locator_a"
        const util::SiPtr< const biol::Atom> end_atom_a( atom_locator_a.LocateAtom( END_STATE));

        // make sure atom could be located in "END_STATE"
        BCL_Assert( end_atom_a.IsDefined(), "end_atom_a is not defined");

        // create atom locator to locate second atom in restraint denoted by "restraint_itr"
        const assemble::LocatorAtom atom_locator_b
        (
          restraint_itr->First().Second().First(), //< chain
          restraint_itr->First().Second().Second(), //< seq id
          restraint_itr->First().Second().Third().GetType() //< atom type
        );

        // get pointer to atom located with "atom_locator_b"
        const util::SiPtr< const biol::Atom> start_atom_b( atom_locator_b.LocateAtom( START_STATE));

        // make sure atom could be located in "START_STATE"
        BCL_Assert( start_atom_b.IsDefined(), "start_atom_b is not defined");

        // get pointer to atom located with "atom_locator_b"
        const util::SiPtr< const biol::Atom> end_atom_b( atom_locator_b.LocateAtom( END_STATE));

        // make sure atom could be located in "END_STATE"
        BCL_Assert( end_atom_b.IsDefined(), "end_atom_b is not defined");

        // get the distance as measured from the starting model
        const double start_distance( linal::Distance( start_atom_a->GetCoordinates(), start_atom_b->GetCoordinates()));

        // get the distance as measured from the ending model
        const double end_distance( linal::Distance( end_atom_a->GetCoordinates(), end_atom_b->GetCoordinates()));

        // calculate the distance change from the start to end model
        const double distance_change( math::Absolute( start_distance - end_distance));

        // message information about the start and end distances and distance change
        BCL_MessageDbg
        (
          "distance change : " +
          util::Format()( restraint_itr->First().First().First())     + " " + //< chain a
          util::Format()( restraint_itr->First().First().Second())    + " " + //< seq id a
          restraint_itr->First().First().Third().GetType().GetName()  + " " + //< atom type a
          "->" +                                                        " " +
          util::Format()( restraint_itr->First().Second().First())    + " " + //< chain b
          util::Format()( restraint_itr->First().Second().Second())   + " " + //< seq id b
          restraint_itr->First().Second().Third().GetType().GetName() +       //< atom type b
          " start_distance " + util::Format()( start_distance) +              //< start distance
          " end_distance " + util::Format()( end_distance) +                  //< end distance
          " distance_change " + util::Format()( distance_change)              //< distance change
        );

        // add the distance change and associated restraint into "distance_changes_and_info"
        distance_changes_and_info.insert
        (
          std::make_pair( distance_change, *restraint_itr)
        );
      }

      // create vector that will hold the restraints ordered by distance change from largest to smallest
      storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > ordered_restraints;

      // iterate through the restraints ordered by distance change from largest to smallest
      // use reverse iterators so the restraints are ordered from largest to smallest
      for
      (
        std::multimap
        <
          double,
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        >::const_reverse_iterator
          itr( distance_changes_and_info.rbegin()), itr_end( distance_changes_and_info.rend());
        itr != itr_end; ++itr
      )
      {
        // add the current restraint to "ordered_restraints"
        ordered_restraints.PushBack( itr->second);
      }

      // end
      return ordered_restraints;
    }

    //! @brief WriteNumberOfDesiredRestraints outputs the desired number of restraints in the desired format to the
    //!         desired filename
    //! @param RESTRAINTS the restraints some of which will be output
    //! @param NUMBER_TO_OUTPUT the number of restraints that will be output
    void RestraintSimulateDistances::WriteNumberOfDesiredRestraints
    (
      const storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > &RESTRAINTS,
      const size_t NUMBER_TO_OUTPUT
    ) const
    {
      // get the number of restraints that will be output
      const size_t num_to_output( RESTRAINTS.GetSize() > NUMBER_TO_OUTPUT ? NUMBER_TO_OUTPUT : RESTRAINTS.GetSize());

      BCL_MessageDbg( "number restraints that will be output : " + util::Format()( num_to_output));

      // get vector with the desired number of restraints
      const storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > restraints_to_write( RESTRAINTS.Begin(), RESTRAINTS.Begin() + num_to_output);

      // create ofstream to write restraints
      io::OFStream write;

      // open "write" to the filename given over the command line
      io::File::MustOpenOFStream( write, m_OutputFile->GetFirstParameter()->GetValue());

      // true if flag was set to write the restraints in rosetta 2.0 format
      if( m_WriteRosettaFormattedRestraints->GetFlag())
      {
        WriteDistanceRestraintsRosettaFormat( restraints_to_write, write);
      }
      // true if flag was set to write the restraints in rosetta 3.0 format
      else if( m_WriteRosettaMiniFormattedRestraints->GetFlag())
      {
        WriteDistanceRestraintsRosettaMiniFormat( restraints_to_write, write);
      }
      // true if flag was set to write restraints in modified CASP RR format
      else if( m_WriteModifiedCASPFormattedRestraints->GetFlag())
      {
        WriteDistanceRestraintsModifiedCASPFormat( restraints_to_write, write);
      }
      else //< write bcl formatted
      {
        restraint::HandlerAtomDistanceAssigned().WriteRestraints( write, restraints_to_write);
      }
    }

    //! @brief GetExplicitlyAllowedRestraints gives the restraint information for restraints based on the residue lists
    //!        passed over the command line. The restraints have to involve residues that are also in the provided list.
    //!        So restraint information is created for all residue pairs between the two command line lists, but then
    //!        restraints are removed that involve residues not contained in the residue list parameter.
    //! @param RESIDUES the list of residues that can be involved in restraints that the restraints are filtered
    //!        according to
    //! @return vector of restraint information
    storage::Vector
    <
      storage::Pair
      <
        storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
      >
    > RestraintSimulateDistances::GetExplicitlyAllowedRestraints
    (
      const util::SiPtrVector< const biol::AABase> &RESIDUES
    ) const
    {
      // create vector which will hold all of the restraint information
      storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > restraint_list;

      // get the filename of the first list of residues to be considered for restraints
      const std::string resi_list_filename_a( m_ExplicitResidueListFilenameA->GetValue());

      // get the filename of the second list of residues to be considered for restraints
      const std::string resi_list_filename_b( m_ExplicitResidueListFilenameB->GetValue());

      // get a collector that will hold the information from "resi_list_filename_a"
      const assemble::CollectorAASpecified collector_a( resi_list_filename_a);

      // get a collector that will hold the information from "resi_list_filename_b"
      const assemble::CollectorAASpecified collector_b( resi_list_filename_b);

      // iterate through the residues in "collector_a"
      for
      (
        storage::List< assemble::LocatorAA>::const_iterator
          resi_list_itr_a( collector_a.GetResidueList().Begin()),
          resi_list_itr_a_end( collector_a.GetResidueList().End());
        resi_list_itr_a != resi_list_itr_a_end;
        ++resi_list_itr_a
      )
      {
        // locate residue in "RESIDUES"
        util::SiPtr< const biol::AABase> resi_a( resi_list_itr_a->Locate( RESIDUES));

        // true if the residue denoted by the locator of "resi_list_itr_a" does not exist in "RESIDUES"
        if( !resi_a.IsDefined())
        {
          // go to next iteration since don't want this residue if it is not in "RESIDUES"
          continue;
        }

        // iterate through the residues in "collector_b"
        for
        (
          storage::List< assemble::LocatorAA>::const_iterator
            resi_list_itr_b( collector_b.GetResidueList().Begin()),
            resi_list_itr_b_end( collector_b.GetResidueList().End());
          resi_list_itr_b != resi_list_itr_b_end;
          ++resi_list_itr_b
        )
        {
          // locate residue in "RESIDUES"
          util::SiPtr< const biol::AABase> resi_b( resi_list_itr_b->Locate( RESIDUES));

          // true if the residue denoted by the locator of "resi_list_itr_b" does not exist in "RESIDUES"
          if( !resi_b.IsDefined())
          {
            // go to next iteration since don't want this residue if it is not in "RESIDUES"
            continue;
          }

          // get the distance between resi_a and resi_b first side chain atoms
          const storage::VectorND< 3, double> distance_upper_lower_bound
          (
            GetDistanceUpperLowerBounds
            (
              resi_a->GetFirstSidechainAtom().GetCoordinates(), resi_b->GetFirstSidechainAtom().GetCoordinates()
            )
          );

          // true if the distance is not defined. For example could be undefined if it is outside the boundaries
          // of distances that should be considered for restraints
          if( !util::IsDefined( distance_upper_lower_bound.First()))
          {
            // continue to next iteration since don't want to add this restraint if distance is undefined
            continue;
          }

          // add the restraint information to "restraint_list"
          restraint_list.PushBack( GetRestraintInformation( *resi_a, *resi_b, distance_upper_lower_bound));
        }
      }

      // return restraint_list
      return restraint_list;
    }

    //! @brief GetExactRestraints gives the restraint information for restraints based on the list of restraints
    //!        passed over the command line. The list contains all pairwise restraints that will be created.
    //!        Restraints are removed that involve residues not contained in the argument parameter residue list
    //!        (which is usually the list of exposed residues).
    //!        The format of the file should be one restraint per line as :
    //!        <chain_col> <resi_seq_id_col> <chain_col> <resi_seq_id_col>
    //!        an example would be :
    //!        'A' 23 'A' 46
    //!        'B' 19 'C' 10
    //! @param RESIDUES the list of residues that can be involved in restraints that the restraints are filtered
    //!        according to
    //! @return vector of restraint information
    storage::Vector
    <
      storage::Pair
      <
        storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
      >
    > RestraintSimulateDistances::GetExactRestraints
    (
      const util::SiPtrVector< const biol::AABase> &RESIDUES
    ) const
    {
      // create vector which will hold all of the restraint information
      storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > restraint_list;

      // get the filename of the file containing the desired restraints
      const std::string restraint_list_filename( m_RestraintListFilename->GetValue());

      // columns for chain and resi ids
      const size_t chain_a_col( m_RestraintListChainAColumn->GetNumericalValue< size_t>());
      const size_t resi_a_col( m_RestraintListResiAColumn->GetNumericalValue< size_t>());
      const size_t chain_b_col( m_RestraintListChainBColumn->GetNumericalValue< size_t>());
      const size_t resi_b_col( m_RestraintListResiBColumn->GetNumericalValue< size_t>());

      io::IFStream read;

      io::File::MustOpenIFStream( read, restraint_list_filename);

      // read in the restraint lines
      storage::Vector< storage::Vector< std::string> > split_lines( util::SplittedStringLineListFromIStream( read));

      // iterate through the restraint lines
      for
      (
        storage::Vector< storage::Vector< std::string> >::const_iterator
          line_itr( split_lines.Begin()), line_itr_end( split_lines.End());
        line_itr != line_itr_end;
        ++line_itr
      )
      {
        // current split line
        const storage::Vector< std::string> &line( *line_itr);

        // skip to next line if this line doesn't have enough columns to access any of the chain or resi ids
        const size_t num_cols( line.GetSize());
        if( num_cols <= chain_a_col || num_cols <= resi_a_col || num_cols <= chain_b_col || num_cols <= resi_b_col)
        {
          BCL_MessageCrt( "skipping line " + util::Format()( line));
          continue;
        }

        // get the chain and resi ids
        std::string chain_a( line( chain_a_col));
        std::string resi_a(  line( resi_a_col));
        std::string chain_b( line( chain_b_col));
        std::string resi_b(  line( resi_b_col));

        // trim the strings
        chain_a = util::TrimString( chain_a);
        resi_a =  util::TrimString( resi_a);
        chain_b = util::TrimString( chain_b);
        resi_b =  util::TrimString( resi_b);

        // if either of the resi ids are not numeric skip to next line or the chains are not one character
        if( !util::IsNumerical( resi_a) || !util::IsNumerical( resi_b) || chain_a.size() != 1 || chain_b.size() != 1)
        {
          BCL_MessageCrt( "skipping line " + util::Format()( line));
          continue;
        }
        // create locators for residues of first restraint
        const assemble::LocatorAA locator_a( chain_a[ 0], util::ConvertStringToNumericalValue< size_t>( resi_a));
        const assemble::LocatorAA locator_b( chain_b[ 0], util::ConvertStringToNumericalValue< size_t>( resi_b));

        // locate resi_a in "RESIDUES"
        util::SiPtr< const biol::AABase> resi_a_ptr( locator_a.Locate( RESIDUES));

        // locate resi_b in "RESIDUES"
        util::SiPtr< const biol::AABase> resi_b_ptr( locator_b.Locate( RESIDUES));

        // true if both residues were found in "RESIDUES"
        if( resi_a_ptr.IsDefined() && resi_b_ptr.IsDefined())
        {
          // get the distance between resi_a and resi_b first side chain atoms
          const storage::VectorND< 3, double> distance_upper_lower_bound
          (
            GetDistanceUpperLowerBounds
            (
              resi_a_ptr->GetFirstSidechainAtom().GetCoordinates(), resi_b_ptr->GetFirstSidechainAtom().GetCoordinates()
            )
          );

          // true if the distance is not defined. For example could be undefined if it is outside the boundaries
          // of distances that should be considered for restraints
          if( util::IsDefined( distance_upper_lower_bound.First()))
          {
            // add the restraint information to "restraint_list"
            restraint_list.PushBack( GetRestraintInformation( *resi_a_ptr, *resi_b_ptr, distance_upper_lower_bound));
          }
        }
      }

      // return restraint_list
      return restraint_list;
    }

    //! @brief GetDistanceUpperLowerBounds gives the distance between the first two coordinates as well as the upper
    //!        and lower bound
    //! @param COORDINATES_A the first coordinate involved in the restraint
    //! @param COORDINATES_B the second coordinate involved in the restraint
    //! @return VectorND< 3> with the distance, upper bound, and lower bound, respectively, of the restraint
    storage::VectorND< 3, double> RestraintSimulateDistances::GetDistanceUpperLowerBounds
    (
      const linal::Vector3D &COORDINATES_A, const linal::Vector3D &COORDINATES_B
    ) const
    {
      // create double "upper_distance_limit" and initialize with "m_UpperLimit"
      const double upper_distance_limit( m_UpperLimit->GetNumericalValue< double>());

      // create double "lower_distance_limit" and initialize with "m_LowerLimit"
      const double lower_distance_limit( m_LowerLimit->GetNumericalValue< double>());

      // create double "distance" and initialize with the distance from "coordinates_a" to "coordinates_b"
      double distance( linal::Distance( COORDINATES_A, COORDINATES_B));

      BCL_Assert( COORDINATES_A.IsDefined(), "coordinates not defined");
      BCL_Assert( COORDINATES_B.IsDefined(), "coordinates not defined");
      BCL_Assert( util::IsDefined( distance), "distance not defined");

      // create the upper bound of the restraint
      const double restraint_upper_bound( distance * 100.0);

      // create the lower bound of the restraint
      double restraint_lower_bound( -100.0);

      // true if "distance" falls outside "upper_distance_limit" and "lower_distance_limit"
      // i.e. true if distance does not meet the distance criteria
      if( distance < lower_distance_limit && distance > upper_distance_limit)
      {
        distance = util::GetUndefinedDouble();
      }

      // true if random distribution amount should be added to distance and the distance is defined, i.e. it falls
      // within the minimum and maximum allowable distances
      if( m_AddDistributionBiasToDistance->GetFlag() && util::IsDefined( distance))
      {
        // get a random value from the distribution
        restraint_lower_bound = GetDistributionBiasAmount();

        // message the bias amount
        BCL_MessageDbg
        (
          "Adding distribution bias amount of " + util::Format()( restraint_lower_bound)
        );

        // add the bias amount to "distance"
        distance += restraint_lower_bound;
      }

      // create vectorND to hold the distance, upper, and lower bounds
      storage::VectorND< 3, double> distance_upper_lower( distance, restraint_upper_bound, restraint_lower_bound);

      // return "distance_upper_lower"
      return distance_upper_lower;
    }

    //! @brief GetRestraintInformation gets the information necessary for a restraint from two residues and the
    //!        distance, upper, and lower bounds
    //! @param RESI_A the first residue involved in the restraint
    //! @param RESI_B the second residue involved in the restraint
    //! @param DISTANCE_UPPER_LOWER_BOUND the distance, upper bound, lower bound, respectively of the restraint
    //! @return the information about the restraint
    storage::Pair
    <
      storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
    > RestraintSimulateDistances::GetRestraintInformation
    (
      const biol::AABase &RESI_A, const biol::AABase &RESI_B,
      const storage::VectorND< 3, double> &DISTANCE_UPPER_LOWER_BOUND
    ) const
    {
      // get the first side chain atom of "RESI_A"
      const biol::Atom &atom_a( RESI_A.GetFirstSidechainAtom());

      // the chain id from resi_a
      const char chain_a( RESI_A.GetChainID());

      // get the seq id from resi_a
      const int seq_id_a( RESI_A.GetSeqID());

      // get the first side chain atom of "RESI_B"
      const biol::Atom &atom_b( RESI_B.GetFirstSidechainAtom());

      // the chain id from resi_b
      const char chain_b( RESI_B.GetChainID());

      // get the seq id from resi_b
      const int seq_id_b( RESI_B.GetSeqID());

      // create pair of restraint information
      const storage::Pair
      <
        storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
      > restraint_info
      (
        storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >
        (
          storage::Triplet< char, int, biol::Atom>( chain_a, seq_id_a, atom_a),
          storage::Triplet< char, int, biol::Atom>( chain_b, seq_id_b, atom_b)
        ),
        DISTANCE_UPPER_LOWER_BOUND
      );

      // get iterator to end of the map
      storage::Map< storage::Pair< char, int>, char>::iterator itr_end( m_ResidueMap.End());

      // create first and second pairs to insert
      storage::Pair< char, int> key_a( restraint_info.First()( 0).First(), restraint_info.First()( 0).Second());
      storage::Pair< char, int> key_b( restraint_info.First()( 1).First(), restraint_info.First()( 1).Second());

      // insert resi_a into map if not already present
      if( m_ResidueMap.Find( key_a) == itr_end)
      {
        m_ResidueMap.InsertElement
        (
          storage::Pair< storage::Pair< char, int>, char>
          (
            key_a,
            RESI_A.GetType()->GetOneLetterCode()
          )
        );
      }

      // insert resi b into map if not already present
      if( m_ResidueMap.Find( key_b) == itr_end)
      {
        m_ResidueMap.InsertElement
        (
          storage::Pair< storage::Pair< char, int>, char>
          (
            key_b,
            RESI_B.GetType()->GetOneLetterCode()
          )
        );
      }

      // return "restraint_info"
      return restraint_info;
    }

    //! @brief creates a shptr vector of Atom Distance Restraints to be written
    //! @param PROTEIN_MODEL where list of residues and restraint distances come from, must be protonated
    //!        and AAComplete
    //! @param NUMBER_DESIRED_RESTRAINTS how many restraints should be generated
    //! @param UPPER_LIMIT_DEVIATION the amount which should be added to the distance found to find the upper bound
    //! @param AA_DISTANCE the smallest sequence distance two aa's can be in sequence for a restraint to be made
    //! @param DISTANCE_RANGE all distances must be in this range to be considered
    //! @return a shptr vector of Atom Distances
    util::ShPtrVector< restraint::AtomDistance> RestraintSimulateDistances::CreateNOERestraints
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const double NUMBER_DESIRED_RESTRAINTS,
      const double UPPER_LIMIT_DEVIATION,
      const size_t AA_DISTANCE,
      const math::Range< double> &DISTANCE_RANGE
    ) const
    {
      // initialize static lower bound for restraints
      static const double s_lower_bound( 1.8);
      const math::Range< double> corrected_range
      (
        std::max( DISTANCE_RANGE.GetMin(), s_lower_bound),
        DISTANCE_RANGE.GetMax()
      );

      // get all atom pairs
      storage::List
      <
        storage::VectorND< 2, storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> > >
      > all_possible_restraints( GetAllProtonPairs( PROTEIN_MODEL, AA_DISTANCE, corrected_range));

      // create size_t "number_restraints" to hold the number of restraint which have been created and initialize to 0
      size_t number_restraints( 0);

      // define a ShPtrVector of Atom Distances to be used for storage
      util::ShPtrVector< restraint::AtomDistance> distance_restraints;

      // continue until the number of restraints has been satisfied or there are no more left to chose from
      while( number_restraints != NUMBER_DESIRED_RESTRAINTS && !all_possible_restraints.IsEmpty())
      {
        // get a random iterator on the list
        storage::List
        <
          storage::VectorND< 2, storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> > >
        >::iterator restraint_itr
        (
          random::GetGlobalRandom().Iterator
          (
            all_possible_restraints.Begin(), all_possible_restraints.End(), all_possible_restraints.GetSize()
          )
        );

        // Get the chain ID
        const char chain_a( restraint_itr->First().Second()->GetChainID());
        const char chain_b( restraint_itr->Second().Second()->GetChainID());

        // Get the Seq ID
        const int seq_id_a( restraint_itr->First().Second()->GetSeqID());
        const int seq_id_b( restraint_itr->Second().Second()->GetSeqID());

        // create SiPtr to atoms
        util::SiPtr< const biol::Atom> &atom_a( restraint_itr->First().First());
        util::SiPtr< const biol::Atom> &atom_b( restraint_itr->Second().First());

        // get the atom types
        const biol::AtomType atom_type_a( atom_a->GetType());
        const biol::AtomType atom_type_b( atom_b->GetType());

        // calculate the distance
        const double calc_distance( linal::Distance( atom_a->GetCoordinates(), atom_b->GetCoordinates()));
        const util::ShPtr< restraint::Distance> distance
        (
          new restraint::Distance( calc_distance, calc_distance + UPPER_LIMIT_DEVIATION, s_lower_bound)
        );

        // create a new restraint
        util::ShPtr< restraint::AtomDistance> sp_restraint
        (
          new restraint::AtomDistance
          (
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
            (
              new restraint::LocatorCoordinatesHydrogen
              (
                chain_a,
                seq_id_a,
                atom_type_a
              )
            ),
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
            (
              new restraint::LocatorCoordinatesHydrogen
              (
                chain_b,
                seq_id_b,
                atom_type_b
              )
            ),
            distance
           )
         );
        distance_restraints.PushBack( sp_restraint);

        ++number_restraints;
        all_possible_restraints.Remove( restraint_itr);
      }
      return distance_restraints;
    }

    //! @brief GetDistributionBiasAmount provides a value based on the distribution provided over the command line
    //! @return double which is the random value obtained from the distribution
    double RestraintSimulateDistances::GetDistributionBiasAmount() const
    {
      // create ShPtr to a Histogram1DDistribution
      static util::ShPtr< random::Histogram1DDistribution> random_distribution;

      // create double which will hold the lower bound of the histogram leading to the biased distribution
      static double lower_bound;

      // create double which will hold the bin size of the histogram leading to the biased distribution
      static double bin_size;

      // true if
      if( !random_distribution.IsDefined())
      {
        // get the filename of the distribution
        const std::string distribution_filename( m_DistributionBiasFilename->GetValue());

        // create histogram to read in the distribution
        math::Histogram distribution;

        // create stream to read from
        io::IFStream read;

        // open "read" to "distribution_filename"
        io::File::MustOpenIFStream( read, distribution_filename);

        // read the histogram from read into "distribution
        read >> distribution;

        // get the lower boundary of "distribution"
        lower_bound = distribution.GetBoundaries().First();

        // get the bin size of "distribution"
        bin_size = distribution.GetBinSize();

        // create a ShPtr to a Histogram1DDistribution created from "distribution"
        const util::ShPtr< random::Histogram1DDistribution> temp_random_distribution
        (
          new random::Histogram1DDistribution( distribution)
        );

        // set "random_distribution" to "temp_random_distribution"
        random_distribution = temp_random_distribution;
      }

      // get a random value from the random distribution
      const double value( random_distribution->DetermineRandomCase( lower_bound, bin_size));

      // return the random value from the distribution
      return value;
    }

    //! @brief returns whether this atom type is valid for ILV-labelling (methyl in ILV or backbone H otherwise)
    //! @param AA_TYPE AA type this atom is associatied with
    //! @param ATOM_TYPE atom type of the atom
    //! @return whether this atom type is valid for ILV-labelling (methyl in ILV or backbone H otherwise)
    bool RestraintSimulateDistances::IsValidILVAtom( const biol::AAType &AA_TYPE, const biol::AtomType &ATOM_TYPE)
    {
      static storage::Map< biol::AAType, storage::Set< biol::AtomType> > s_methyl_map;
      if( s_methyl_map.IsEmpty())
      {
        s_methyl_map[ biol::GetAATypes().ILE].Insert( biol::GetAtomTypes().HD11);
        s_methyl_map[ biol::GetAATypes().ILE].Insert( biol::GetAtomTypes().HD12);
        s_methyl_map[ biol::GetAATypes().ILE].Insert( biol::GetAtomTypes().HD13);

        s_methyl_map[ biol::GetAATypes().LEU].Insert( biol::GetAtomTypes().HD11);
        s_methyl_map[ biol::GetAATypes().LEU].Insert( biol::GetAtomTypes().HD12);
        s_methyl_map[ biol::GetAATypes().LEU].Insert( biol::GetAtomTypes().HD13);
        s_methyl_map[ biol::GetAATypes().LEU].Insert( biol::GetAtomTypes().HD21);
        s_methyl_map[ biol::GetAATypes().LEU].Insert( biol::GetAtomTypes().HD22);
        s_methyl_map[ biol::GetAATypes().LEU].Insert( biol::GetAtomTypes().HD23);

        s_methyl_map[ biol::GetAATypes().VAL].Insert( biol::GetAtomTypes().HG11);
        s_methyl_map[ biol::GetAATypes().VAL].Insert( biol::GetAtomTypes().HG12);
        s_methyl_map[ biol::GetAATypes().VAL].Insert( biol::GetAtomTypes().HG13);
        s_methyl_map[ biol::GetAATypes().VAL].Insert( biol::GetAtomTypes().HG21);
        s_methyl_map[ biol::GetAATypes().VAL].Insert( biol::GetAtomTypes().HG22);
        s_methyl_map[ biol::GetAATypes().VAL].Insert( biol::GetAtomTypes().HG23);
      }

      // if is ILV type
      if( s_methyl_map.Has( AA_TYPE))
      {
        // return true if it is a methyl
        if( s_methyl_map[ AA_TYPE].Contains( ATOM_TYPE))
        {
          return true;
        }
      }

      // true if backbone H
      return ATOM_TYPE == biol::GetAtomTypes().H;
    }

    //! @brief takes a protein model and returns all the atom pairs that have distances within the range
    //! @param PROTEIN_MODEL where list of residues and restraint distances come from, must be protonated
    //!        and AAComplete
    //! @param AA_DISTANCE the smallest sequence distance two aa's can be in sequence for a restraint to be made
    //! @param DISTANCE_RANGE all distances must be in this range to be considered
    //! @return all the atom pairs that have distances within the range
    storage::List
    <
      storage::VectorND< 2, storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> > >
    > RestraintSimulateDistances::GetAllProtonPairs
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const size_t AA_DISTANCE,
      const math::Range< double> &DISTANCE_RANGE
    ) const
    {
      // initialize return list
      storage::List
      <
        storage::VectorND< 2, storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> > >
      > all_possible_restraints;

      // create SiPtrVector "amino_acids" and initialize with the amino acids in the SSEs of "protein_model"
      const util::SiPtrVector< const biol::AABase> amino_acids( PROTEIN_MODEL.GetAminoAcids());

      // create set of protons to exclude (i.e. N-term and C-term)
      static const storage::Set< biol::AtomType> s_excluded_protons
      (
        biol::GetAtomTypes().H1,
        biol::GetAtomTypes().H2,
        biol::GetAtomTypes().H3,
        biol::GetAtomTypes().HXT
      );

      const bool ilv_only( m_SimulateNMRDistanceRestraint->GetParameterList()( 1)->GetWasSetInCommandLine());

      // iterate through "amino_acids" to build up "all_possible_restraints"
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          itr_a( amino_acids.Begin()), itr_end( amino_acids.End());
        itr_a != itr_end;
        ++itr_a
      )
      {
        // iterate through "amino_acids" to build up "all_possible_restraints"
        for
        (
          util::SiPtrVector< const biol::AABase>::const_iterator itr_b( itr_a + 1);
          itr_b != itr_end;
          ++itr_b
        )
        {
          // if the sequence distance is large enough (or on different chains)
          if
          (
            math::Absolute( ( *itr_a)->GetSeqID() - ( *itr_b)->GetSeqID()) >= int( AA_DISTANCE) ||
              ( *itr_a)->GetChainID() != ( *itr_b)->GetChainID()
          )
          {
            // get the hydrogen atoms
            const util::SiPtrVector< const biol::Atom> h_atoms_a
            (
              ( *itr_a)->GetAtoms( chemistry::GetElementTypes().e_Hydrogen)
            );
            const util::SiPtrVector< const biol::Atom> h_atoms_b
            (
              ( *itr_b)->GetAtoms( chemistry::GetElementTypes().e_Hydrogen)
            );

            // iterate through the h atoms
            for
            (
              util::SiPtrVector< const biol::Atom>::const_iterator atom_itr_a( h_atoms_a.Begin()),
                atom_itr_a_end( h_atoms_a.End());
              atom_itr_a != atom_itr_a_end; ++atom_itr_a
            )
            {
              // if only considering ILV atoms
              if( ilv_only)
              {
                // if atom is not ILV or backbone
                if( !IsValidILVAtom( ( *itr_a)->GetType(), ( *atom_itr_a)->GetType()))
                {
                  // skip this atom
                  continue;
                }
              }

              for
              (
                util::SiPtrVector< const biol::Atom>::const_iterator atom_itr_b( h_atoms_b.Begin()),
                  atom_itr_b_end( h_atoms_b.End());
                atom_itr_b != atom_itr_b_end; ++atom_itr_b
              )
              {
                // if only considering ILV atoms
                if( ilv_only)
                {
                  // if atom is not ILV or backbone
                  if( !IsValidILVAtom( ( *itr_b)->GetType(), ( *atom_itr_b)->GetType()))
                  {
                    // skip this atom
                    continue;
                  }
                }

                // get the distance
                const double distance( linal::Distance( ( *atom_itr_a)->GetCoordinates(), ( *atom_itr_b)->GetCoordinates()));

                // if the euclidean distance is within the range
                if
                (
                  util::IsDefined( distance) && DISTANCE_RANGE.IsWithin( distance) &&
                  !s_excluded_protons.Contains( ( *atom_itr_a)->GetType()) &&
                  !s_excluded_protons.Contains( ( *atom_itr_b)->GetType())
                )
                {
                  // add the atoms to the list
                  all_possible_restraints.PushBack
                  (
                    storage::VectorND
                    <
                      2,
                      storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
                    >
                    (
                      storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
                      (
                        *atom_itr_a,
                        *itr_a
                      ),
                      storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
                      (
                        *atom_itr_b,
                        *itr_b
                      )
                    )
                  );
                } // distance in range
              } // atom b for loop
            } // atom a for loop
          } // sequence distance is large enough
        } // aa b for loop
      } // aa a for loop

      // end
      return all_possible_restraints;
    }

    // default constructor
    RestraintSimulateDistances::RestraintSimulateDistances() :
      m_PDBFilename
      (
        new command::FlagStatic
        (
          pdb::GetDefaultFileExtension(),
          "\tpdb file for which restraints will be created",
          command::Parameter
          (
            "pdb_filename",
            "\tfilename for which restraints will be calculated",
              command::ParameterCheckExtension( ".pdb"),
            ""
          )
        )
      ),
      m_SimulateDistanceRestraint
      (
        new command::FlagStatic
        (
          "simulate_distance_restraints",
          "\tIndicates that distance restraints should be simulated for the given pdb file"
        )
      ),
      m_SkipUndefinedAAs
      (
        new command::FlagStatic
        (
          "skip_undefined_aas",
          "\tIndicates that simulated distance restraints will not consider undefined AAs in the given pdb file"
        )
      ),
      m_SimulateNMRDistanceRestraint
      (
        new command::FlagStatic
        (
          "simulate_nmr_distance_restraints",
          "\tIndicates that nmr distance restraints should be simulated for the given pdb file, which must be "
          "protonated and be read in as AAComplete",
          util::ShPtrVector< command::ParameterInterface>::Create
          (
            util::CloneToShPtr
            (
              command::Parameter
              (
                "aa_distance_apart",
                "the distance apart in sequence two AA must be in order to make a restraint",
                "5"
              )
            ),
            util::CloneToShPtr
            (
              command::Parameter
              (
                "ILV",
                "only use ILE, LEU, and VAL for generating NOE restraints",
                ""
              )
            )
          )
        )
      ),
      m_NumberDesiredRestraints
      (
        new command::FlagStatic
        (
          "number_of_restraints", "Flag for the specifiying an exact number of desired restraints.",
           command::Parameter
          (
            "number_of_desired_restraints",
            "The number of restraints which are desired to be gotten.",
            util::Format()( std::numeric_limits< size_t>::max())
          )
        )
      ),
      m_Limit
      (
        new command::FlagStatic
        (
          "distance_limits",
          "The limits that a measured distance should have."
        )
      ),
      m_UpperLimit
      (
        new command::Parameter
        (
          "upper_limit",
          "\tThe maximum measured distance that should be used to create a restraint.",
          util::Format()( std::numeric_limits< double>::max())
        )
      ),
      m_LowerLimit
      (
        new command::Parameter
        (
          "lower_limit",
          "\tThe minimum measured distance that should be used to create a restraint.",
          "0"
        )
      ),
      m_MinimumSeparation
      (
        new command::FlagStatic
        (
          "minimum_separation",
          "minimum sequence separation to be used for restraint selection, default 0",
          command::Parameter
          (
            "minimum_separation",
            "minimum sequence separation to be used for restraint selection, default 0",
            "0"
          )
        )
      ),
      m_NeighborCountLimit
      (
        new command::FlagStatic
        (
          "neighbor_count_limit",
          "The maximum neighbor count that a residue can have and still be considered for a restraint.",
           command::Parameter
          (
            "neighbor_count",
            "The maximum neighbor count that a residue can have and still be a part of a restraint.",
            util::Format()( std::numeric_limits< double>::max())
          )
        )
      ),
      m_OutputFile
      (
        new command::FlagStatic
        (
          "output_file",
          "Path and name of the output file which will hold the simulated restraints.",
          command::Parameter
          (
            "string",
            "\tpath and name of the outputfile",
            "simulated.cst_dist"
          )
        )
      ),
      m_ConsiderExplicitResidues
      (
        new command::FlagStatic
        (
          "consider_residues",
          "only residues between specified lists should be considered for restraints"
        )
      ),
      m_ExplicitResidueListFilenameA
      (
        new command::Parameter
        (
          "list_filename_a",
          "\tThe first list of residues. Format \"<chain_id> <seq_id>\" one per line",
          "residues_a.ls"
        )
      ),
      m_ExplicitResidueListFilenameB
      (
        new command::Parameter
        (
          "list_filename_b",
          "\tThe second list of residues. Format \"<chain_id> <seq_id>\" one per line",
          "residues_b.ls"
        )
      ),
      m_OrderByDistanceChange
      (
        new command::FlagStatic
        (
          "order_by_distance_change",
          "restraints involving residues changing the most between two states will be given"
        )
      ),
      m_EndStatePDBFilename
      (
        new command::Parameter
        (
          "state_b_pdb_filename",
          "\tthe pdb with the second state",
          "second_state.pdb"
        )
      ),
      m_WriteRosettaFormattedRestraints
      (
        new command::FlagStatic
        (
          "write_rosetta_restraints",
          "restraints will be written in rosetta 2.0 format"
        )
      ),
      m_AddDistributionBiasToDistance
      (
        new command::FlagStatic
        (
          "add_distance_uncertainty",
          "a random amount should be added to the first side chain atom distance between two residues. This can be used to simulate the uncertainty in a epr distance measurement."
        )
      ),
      m_DistributionBiasFilename
      (
        new command::Parameter
        (
          "histogram_distribution_filename",
          "\tthe file which contains the bcl formatted histogram which will be used as the distribution",
          "distribution.bcl_histogram"
        )
      ),
      m_WriteRosettaMiniFormattedRestraints
      (
        new command::FlagStatic
        (
          "write_rosetta_mini_restraints",
          "restraints will be written in rosetta 3.0 format"
        )
      ),
      m_WriteModifiedCASPFormattedRestraints
      (
        new command::FlagStatic
        (
          "write_modified_CASP_restraints",
          "restraints will be written in modified CASP format"
        )
      ),
      m_RestraintListProvided
      (
        new command::FlagStatic
        (
          "restraint_list",
          "a list of desired restraints will be provided"
        )
      ),
      m_RestraintListFilename
      (
        new command::Parameter
        (
          "restraint_list_filename",
          "\tthe file which contains the list of restraints. Restraints per line formatted as : <chain_col> <resi_seq_id_col> <chain_col> <resi_seq_id_col>",
          "restraints.ls"
        )
      ),
      m_RestraintListChainAColumn
      (
        new command::Parameter
        (
          "chain_a_column",
          "\tcolumn of the file which specifies the chain of the first residue. Columns starts at 0",
          "0"
        )
      ),
      m_RestraintListResiAColumn
      (
        new command::Parameter
        (
          "resi_a_column",
          "\tcolumn of the file which specifies the residue number of the first residue. Columns starts at 0",
          "1"
        )
      ),
      m_RestraintListChainBColumn
      (
        new command::Parameter
        (
          "chain_b_column",
          "\tcolumn of the file which specifies the chain of the second residue. Columns starts at 0",
          "2"
        )
      ),
      m_RestraintListResiBColumn
      (
        new command::Parameter
        (
          "resi_b_column",
          "\tcolumn of the file which specifies the residue number of the second residue. Columns starts at 0",
          "3"
        )
      ),
      m_FractionDesiredRestraints
      (
        new command::FlagStatic
        (
          "num_restraint_fraction",
          "flag specifying that the number of restraints should be a fraction of the residues in the protein"
        )
      ),
      m_FractionDesiredRestraintsNumber
      (
        new command::Parameter
        (
          "fraction",
          "\tThe fraction of restraints compared to the number of residues. e.g. 0.1 = 10 restraints for 100 residue protein",
          "0.1"
        )
      )
    {
      // attach parameters to flags
      m_Limit->PushBack( m_UpperLimit);
      m_Limit->PushBack( m_LowerLimit);
      m_ConsiderExplicitResidues->PushBack( m_ExplicitResidueListFilenameA);
      m_ConsiderExplicitResidues->PushBack( m_ExplicitResidueListFilenameB);
      m_OrderByDistanceChange->PushBack( m_EndStatePDBFilename);
      m_AddDistributionBiasToDistance->PushBack( m_DistributionBiasFilename);
      m_RestraintListProvided->PushBack( m_RestraintListFilename); //< filename containing restraints
      m_RestraintListProvided->PushBack( m_RestraintListChainAColumn);
      m_RestraintListProvided->PushBack( m_RestraintListResiAColumn);
      m_RestraintListProvided->PushBack( m_RestraintListChainBColumn);
      m_RestraintListProvided->PushBack( m_RestraintListResiBColumn);
    }

    const ApplicationType RestraintSimulateDistances::RestraintSimulateDistances_Instance
    (
      GetAppGroups().AddAppToGroup( new RestraintSimulateDistances(), GetAppGroups().e_Restraint)
    );

  } // namespace app
} // namespace bcl
