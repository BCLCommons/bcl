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
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "nmr/bcl_nmr_rosetta_noe_handler.h"
#include "nmr/bcl_nmr_rosetta_rdc_handler.h"
#include "nmr/bcl_nmr_star_noe_handler.h"
#include "nmr/bcl_nmr_star_rdc_handler.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"
#include "score/bcl_score_restraint_nmr_distance_interface.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NmrFileConvert
    //! @brief app is for switching between ROSETTA and star nmr restraint files.
    //!        1. Can also specify number of restraints to be created, which will get a random subset from the full set
    //!           of the restraints.
    //!        2. The sequence offset function allows the user to compensate for inconsistencies in sequence ID
    //!           numbering due to PDB and restraint files that have numbering starting at some number other than 1.  A
    //!           PDB should only be passed to the same set of restraints one time to fix this problem.
    //!        3. Can also convert all atoms to back bone atoms which will change all non backbone atoms to CB and add
    //!           one angstrom for each bond back to the backbone atom to the distance.  This is necessary for ROSETTA
    //!           Abinitio folding.
    //!
    //! @author akinlr
    //! @date 06/28/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API NmrFileConvert :
      public Interface
    {
    private:

    //////////
    // data //
    //////////

      //! command line parameters
      util::ShPtr< command::FlagStatic> m_InputFileFlag;
      util::ShPtr< command::Parameter> m_InputFileParam;
      util::ShPtr< command::Parameter> m_InputFileTypeParam;
      util::ShPtr< command::Parameter> m_InputDataTypeParam;
      util::ShPtr< command::FlagStatic> m_OutputFileFlag;
      util::ShPtr< command::Parameter> m_OutputFileParam;
      util::ShPtr< command::Parameter> m_OutputFileTypeParam;
      util::ShPtr< command::FlagStatic> m_PdbFileFlag;
      util::ShPtr< command::FlagStatic> m_RandomSubsetFlag;
      util::ShPtr< command::FlagStatic> m_SequenceSeparationFlag;
      util::ShPtr< command::FlagStatic> m_GetCentroidDistancesFlag;

      //! flag for KB histogram prefix
      util::ShPtr< command::FlagStatic> m_KBPotentialsPrefixFlag;

      //! flag for rosetta score weight
      util::ShPtr< command::FlagStatic> m_RosettaWeightFlag;

      //! flag for whether to adjust sequence
      util::ShPtr< command::FlagInterface> m_AdjustSeqIDFlag;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      NmrFileConvert();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      NmrFileConvert *Clone() const
      {
        return new NmrFileConvert( *this);
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

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief GenerateProteinModel creates and returns a protein model based on a command flag
      //! @param PDB_STRING is the string which gives the pdb to create the protein model from
      //! @return ProteinModel which was created from the pdb given by "PDB_STRING"
      assemble::ProteinModel GenerateProteinModel( const std::string &PDB_STRING) const;

    private:

      //! @brief modify distances and atom names based on distance from the backbone
      //! @return modified NOEs
      util::ShPtrVector< restraint::AtomDistance> ConvertToCentroid
      (
        const util::ShPtrVector< restraint::AtomDistance> &ORIGINAL_RESTRAINTS
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

    private:

      static const ApplicationType NmrFileConvert_Instance;

    }; //class NmrFileConvert

  /////////////////
  // data access //
  /////////////////

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> NmrFileConvert::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add flags
      sp_cmd->AddFlag( m_InputFileFlag);
      sp_cmd->AddFlag( m_OutputFileFlag);
      sp_cmd->AddFlag( m_PdbFileFlag);
      sp_cmd->AddFlag( m_RandomSubsetFlag);
      sp_cmd->AddFlag( m_SequenceSeparationFlag);
      sp_cmd->AddFlag( m_GetCentroidDistancesFlag);
      sp_cmd->AddFlag( m_KBPotentialsPrefixFlag);
      sp_cmd->AddFlag( m_RosettaWeightFlag);
      sp_cmd->AddFlag( m_AdjustSeqIDFlag);

      // add pdb flags
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    int NmrFileConvert::Main() const
    {
      // create a stream to read in files
      io::IFStream read;

      // store all of the flags
      const std::string input_file( m_InputFileParam->GetValue());
      const std::string input_file_type( m_InputFileTypeParam->GetValue());
      const std::string data_type( m_InputDataTypeParam->GetValue());
      const std::string output_file_type( m_OutputFileTypeParam->GetValue());
      const std::string output_file( m_OutputFileParam->GetValue());
      size_t subset_nr( m_RandomSubsetFlag->GetFirstParameter()->GetNumericalValue< size_t>());
      const size_t sequence_sep( m_SequenceSeparationFlag->GetFirstParameter()->GetNumericalValue< size_t>());

      // initialize the protein model so it can be used if a pdb file is passed to adjust the numbering in the
      // restraint file
      assemble::ProteinModel protein_model;
      size_t seq_offset( 0);

      if( m_PdbFileFlag->GetFlag())
      {
        // create ProteinModel "protein_model" and initialize with protein of "pdb_filename"
        protein_model = GenerateProteinModel( m_PdbFileFlag->GetFirstParameter()->GetValue());

        const util::SiPtrVector< const biol::AABase> aas( protein_model.GetAminoAcids());

        // if the seq offset should be set
        if( m_AdjustSeqIDFlag->GetFlag())
        {
          // get the sequence offset value from the pdb id of the first AA
          seq_offset = aas.FirstElement()->GetPdbID() - aas.FirstElement()->GetSeqID();

          BCL_MessageStd( "Applying sequence offset: " + util::Format()( seq_offset));
        }

        // if the subset number is undefined
        if( !util::IsDefined( subset_nr))
        {
          // set it to the number of residues in SSEs
          subset_nr = aas.GetSize();
        }
      }

      // only make the type of restraint necessary based on the flag output
      if( data_type == "noe")
      {
        // initialize the Handlers
        nmr::RosettaNOEHandler rosetta_handler
        (
          sequence_sep,
          seq_offset,
          m_KBPotentialsPrefixFlag->GetFirstParameter()->GetValue(),
          m_RosettaWeightFlag->GetFirstParameter()->GetNumericalValue< double>()
        );
        nmr::StarNOEHandler star_handler( "", sequence_sep, seq_offset);

        // read in the file
        io::File::MustOpenIFStream( read, input_file);

        // initialize an empty vector of noes
        util::ShPtrVector< restraint::AtomDistance> noes;

        // determine how to read the file based on file and data types
        if( input_file_type == "star")
        {
          // read in using the star handler if it's the star file type
          noes = star_handler.ReadRestraints( read);
        }
        else
        {
          // read in using the rosetta handler if it's the rosetta file type
          noes = rosetta_handler.ReadRestraints( read);
        }

        // close the file
        io::File::CloseClearFStream( read);

        BCL_MessageStd( "number of NOEs originally is: " + util::Format()( noes.GetSize()));

        // initialize another empty vector to store the cleaned up noes
        util::ShPtrVector< restraint::AtomDistance> clean_noes;

        // if a protein model was given
        if( m_PdbFileFlag->GetFlag())
        {
          // iterate over the restraints
          for
          (
            util::ShPtrVector< restraint::AtomDistance>::const_iterator itr( noes.Begin()), end_itr( noes.End());
            itr != end_itr; ++itr
          )
          {
            // locate the two atom types
            const biol::AtomType type_a( ( *itr)->GetData().First()->LocateAtomCopy( protein_model).GetType());
            const biol::AtomType type_b( ( *itr)->GetData().Second()->LocateAtomCopy( protein_model).GetType());

            // if the restraints can be located in the protein model
            if( type_a.IsDefined() && type_b.IsDefined())
            {
              // add the restraint
              clean_noes.PushBack
              (
                util::ShPtr< restraint::AtomDistance>
                (
                  new restraint::AtomDistance
                  (
                    util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
                    (
                      new restraint::LocatorCoordinatesHydrogen
                      (
                        ( *itr)->GetData().First()->GetChainID(),
                        ( *itr)->GetData().First()->GetSeqID(),
                        type_a.GetName()
                      )
                    ),
                    util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
                    (
                      new restraint::LocatorCoordinatesHydrogen
                      (
                        ( *itr)->GetData().Second()->GetChainID(),
                        ( *itr)->GetData().Second()->GetSeqID(),
                        type_b.GetName()
                      )
                    ),
                    ( *itr)->GetDistance()
                  )
                )
              );
            }
            // restraint not found in model
            else
            {
              BCL_MessageStd
              (
                "Following restraint not found in protein model: " + ( *itr)->GetIdentification()
              );
            }
          }

          // store the clean noes
          noes = clean_noes;
        }

        BCL_MessageStd( "number of NOEs remaining is: " + util::Format()( noes.GetSize()))

        // if Get Centroid flag is given
        if( m_GetCentroidDistancesFlag->GetFlag())
        {
          // convert the NOE set to Centroid mode
          noes = ConvertToCentroid( noes);

          BCL_MessageStd( "number of NOEs remaining after centroid mode is: " + util::Format()( noes.GetSize()))
        }

        // if flag for getting random subset is given
        if( m_RandomSubsetFlag->GetFlag())
        {
          // store the number of restraints that have been stored in the newest file
          size_t curr_nr_restraints( 0);

          // initialize a new vector of restraints to be retrieved
          util::ShPtrVector< restraint::AtomDistance> random_subset_noes;

          // store restraints until the number in the subset is fulfilled or there are no noes left
          while( ( curr_nr_restraints < subset_nr) && ( !noes.IsEmpty()))
          {
            // take an element from the original subset and store it in the new subset
            random_subset_noes.PushBack( noes.RemoveRandomElement());

            curr_nr_restraints = random_subset_noes.GetSize();
          }

          // store the random subset as the set to be returned
          noes = random_subset_noes;
        }

        // prepare to write the new file
        io::OFStream write;
        io::File::MustOpenOFStream( write, output_file);

        // store the file based on the file type given
        if( output_file_type == "star")
        {
          // store as a star file
          star_handler.WriteRestraints( write, noes);
        }
        else
        {
          // if a pdb file was given
          if( m_PdbFileFlag->GetFlag())
          {
            // iterate over the restraints
            for
            (
              util::ShPtrVector< restraint::AtomDistance>::iterator itr( noes.Begin()), itr_end( noes.End());
              itr != itr_end; ++itr
            )
            {
              // get the locators
              const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> sp_locator_a( ( *itr)->GetData().First());
              const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> sp_locator_b( ( *itr)->GetData().Second());

              // get the pdb atom types
              std::string pdb_type_a;
              std::string pdb_type_b;
              if( m_KBPotentialsPrefixFlag->GetFlag())
              {
                // locate the two atom types
                pdb_type_a = ( *itr)->GetData().First()->LocateAtomCopy( protein_model).GetType().GetName();
                pdb_type_b = ( *itr)->GetData().Second()->LocateAtomCopy( protein_model).GetType().GetName();
              }
              else
              {
                pdb_type_a = sp_locator_a->LocateAA( protein_model)->GetType()->GetPDBAtomName
                (
                  sp_locator_a->LocateAtomCopy( protein_model).GetType()
                );
                pdb_type_b =
                sp_locator_b->LocateAA( protein_model)->GetType()->GetPDBAtomName
                (
                  sp_locator_b->LocateAtomCopy( protein_model).GetType()
                );
              }

              // update the noe
              *itr = util::ShPtr< restraint::AtomDistance>
              (
                new restraint::AtomDistance
                (
                  util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
                  (
                    new restraint::LocatorCoordinatesHydrogen
                    (
                      sp_locator_a->GetChainID(),
                      sp_locator_a->GetSeqID(),
                      pdb_type_a
                    )
                  ),
                  util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
                  (
                    new restraint::LocatorCoordinatesHydrogen
                    (
                      sp_locator_b->GetChainID(),
                      sp_locator_b->GetSeqID(),
                      pdb_type_b
                    )
                  ),
                  ( *itr)->GetDistance()
                )
              );
            }
          }

          // store as a Rosetta file
          rosetta_handler.WriteRestraints( write, noes);
        }
        io::File::CloseClearFStream( write);
      }

      // RDC DATA

      // only make the type of restraint necessary based on the flag output
      if( data_type == "rdc")
      {
        // initialize the Handlers
        nmr::RosettaRDCHandler rosetta_handler( seq_offset);
        nmr::StarRDCHandler star_handler( "", seq_offset);

        // open file to read it in
        io::File::MustOpenIFStream( read, input_file);

        // initialize and empty vector of rdcs
        restraint::RDC rdcs;

        // decide how to read the file based on file and data types
        if( input_file_type == "star")
        {
          rdcs = star_handler.ReadRestraints( read);
        }
        else
        {
          rdcs = rosetta_handler.ReadRestraints( read);
        }

        // close the file once it has been read in
        io::File::CloseClearFStream( read);

        // if a protein model was given
        if( m_PdbFileFlag->GetFlag())
        {
          // initialize cleaned rdcs
          restraint::RDC clean_rdcs;

          // iterate over the restraints
          for
          (
            storage::Vector< storage::Triplet< restraint::DataPairwise, double, double> >::const_iterator
              itr( rdcs.GetData().Begin()), itr_end( rdcs.GetData().End());
            itr != itr_end; ++itr
          )
          {
            // locate the two atom types
            const biol::AtomType type_a( itr->First().First()->LocateAtomCopy( protein_model).GetType());
            const biol::AtomType type_b( itr->First().Second()->LocateAtomCopy( protein_model).GetType());

            // if the restraints can be located in the protein model
            if( type_a.IsDefined() && type_b.IsDefined())
            {
              // add the new rdc data w/ the pdb atom type
              clean_rdcs.PushBack
              (
                restraint::LocatorCoordinatesHydrogen
                (
                  itr->First().First()->GetChainID(),
                  itr->First().First()->GetSeqID(),
                  type_a
                ),
                restraint::LocatorCoordinatesHydrogen
                (
                  itr->First().Second()->GetChainID(),
                  itr->First().Second()->GetSeqID(),
                  type_b
                ),
                itr->Second(),
                itr->Third()
              );
            }
            // restraint not found in model
            else
            {
              BCL_MessageStd
              (
                "Following restraint not found in protein model: " + itr->First().First()->GetIdentification() +
                " : " + itr->First().Second()->GetIdentification()
              );
            }
          }

          // store the clean noes
          rdcs = clean_rdcs;
        }

        // take a random subset of the restraints if flag is given
        if( m_RandomSubsetFlag->GetFlag())
        {
          // initialize a value for the current new subset
          size_t curr_nr_restraints( 0);

          // store the number of RDCs
          size_t restraint_nr( rdcs.GetData().GetSize());

          // initialize a new set for the random subset
          restraint::RDC random_subset_rdcs;

          // get the data from the original rdcs so random parts of it can be removed
          storage::Vector< storage::Triplet< restraint::DataPairwise, double, double> > temp_rdcs( rdcs.GetData());

          // store RDCs until the random subset number has been reached or all of the restraints have been stored
          while( curr_nr_restraints <= subset_nr && curr_nr_restraints <= restraint_nr)
          {
            // initialize a random
            storage::Triplet< restraint::DataPairwise, double, double> random_rdc = temp_rdcs.RemoveRandomElement();

            // store the random rdcs
            random_subset_rdcs.PushBack
            (
              *( random_rdc.First().First()),
              *( random_rdc.First().Second()),
              random_rdc.Second(),
              random_rdc.Third()
            );

            // increase the number of restraints
            curr_nr_restraints++;
          }

          // store the random subset of rdcs as the set of rdcs
          rdcs = random_subset_rdcs;
        }

        // prepare to write the new file
        io::OFStream write;
        io::File::MustOpenOFStream( write, output_file);

        // write the rdc file based on whether or not star or Rosetta is specified
        if( output_file_type == "star")
        {
          star_handler.WriteRestraints( write, rdcs);
        }
        else
        {
          // if a pdb file was given
          if( m_PdbFileFlag->GetFlag())
          {
            // create new vector
            restraint::RDC pdb_rdcs;

            // iterate over the restraints
            for
            (
              storage::Vector< storage::Triplet< restraint::DataPairwise, double, double> >::const_iterator
                itr( rdcs.GetData().Begin()), itr_end( rdcs.GetData().End());
              itr != itr_end; ++itr
            )
            {
              // dynamic cast the locators to LocatorCoordinatesHydrogen so that the string can be accessed
              const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> sp_locator_a( itr->First().First());
              const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> sp_locator_b( itr->First().Second());

              // get the pdb atom types
              const std::string pdb_type_a
              (
                sp_locator_a->LocateAA( protein_model)->GetType()->GetPDBAtomName
                (
                  sp_locator_a->LocateAtomCopy( protein_model).GetType()
                )
              );
              const std::string pdb_type_b
              (
                sp_locator_b->LocateAA( protein_model)->GetType()->GetPDBAtomName
                (
                  sp_locator_b->LocateAtomCopy( protein_model).GetType()
                )
              );

              // add the new rdc data w/ the pdb atom type
              pdb_rdcs.PushBack
              (
                restraint::LocatorCoordinatesHydrogen
                (
                  sp_locator_a->GetChainID(),
                  sp_locator_a->GetSeqID(),
                  pdb_type_a
                ),
                restraint::LocatorCoordinatesHydrogen
                (
                  sp_locator_b->GetChainID(),
                  sp_locator_b->GetSeqID(),
                  pdb_type_b
                ),
                itr->Second(),
                itr->Third()
              );
            }

            // store the rdcs
            rdcs = pdb_rdcs;
          }

          // write the rdcs
          rosetta_handler.WriteRestraints( write, rdcs);
        }
        io::File::CloseClearFStream( write);
      }

      return 0;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief GenerateProteinModel creates and returns a protein model based on a command flag
    //! @param PDB_STRING is the string which gives the pdb to create the protein model from
    //! @return ProteinModel which was created from the pdb given by "PDB_STRING"
    assemble::ProteinModel NmrFileConvert::GenerateProteinModel( const std::string &PDB_STRING) const
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

    //! @brief modify distances and atom names based on distance from the backbone
    //! @return modified NOEs
    util::ShPtrVector< restraint::AtomDistance>
      NmrFileConvert::ConvertToCentroid( const util::ShPtrVector< restraint::AtomDistance> &ORIGINAL_RESTRAINTS) const
    {
      // initialize a new vector for storing the modified restraints
      util::ShPtrVector< restraint::AtomDistance> new_restraints;

      // iterate through the old list to modify each atom according to its location from the backbone atom
      for
      (
        util::ShPtrVector< restraint::AtomDistance>::const_iterator itr( ORIGINAL_RESTRAINTS.Begin()),
          itr_end( ORIGINAL_RESTRAINTS.End());
        itr != itr_end; itr++
      )
      {
        // dynamic cast the locators to LocatorCoordinatesHydrogen so that the string can be accessed
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_a( ( *itr)->GetData().First());
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_b( ( *itr)->GetData().Second());

        // store the current atom types
        biol::AtomType atom_a( biol::GetAtomTypes().TypeFromPDBAtomName( sp_locator_a->GetAtomTypeString()));
        biol::AtomType atom_b( biol::GetAtomTypes().TypeFromPDBAtomName( sp_locator_b->GetAtomTypeString()));

        // Get Atom Distance from CB before the atom is changed to CB (values should come out 0 if they are bacbone atoms
        const size_t atom_a_dist( score::RestraintNMRDistanceInterface::GetBondsFromCB( atom_a));
        const size_t atom_b_dist( score::RestraintNMRDistanceInterface::GetBondsFromCB( atom_b));

        // Get the Distance so we know how to modify it
        const double final_dist( ( *itr)->GetDistance()->GetDistance() + double( atom_a_dist) + double( atom_b_dist));

        // Tell user what the modified distance is
        BCL_MessageCrt
        (
          "NOE distance has been changed to " + util::Format()( final_dist) +
            " from " + util::Format()( ( *itr)->GetDistance()->GetDistance()) + " to compensate for Centroid mode"
        )

        // modify the upper distance to make it a reasonable value in comparison to the modified distance
        const double final_upper( ( *itr)->GetDistance()->UpperBound() + double( atom_a_dist) + double( atom_b_dist));

        // store the distances
        const util::ShPtr< restraint::Distance> distance
        (
          new restraint::Distance( final_dist, final_upper, ( *itr)->GetDistance()->LowerBound())
        );

        // get all of the backbone atom types to compare with the current atom type
        bool a_is_backbone
        (
          atom_a.GetName() == "HA" ||
          atom_a.GetName() == "H" ||
          atom_a.GetName() == "CA" ||
          atom_a.GetName() == "CB" ||
          atom_a.GetName() == "N" ||
          atom_a.GetName() == "C" ||
          atom_a.GetName() == "HA2" ||
          atom_a.GetName() == "HA3" ||
          atom_a == biol::GetAtomTypes().e_Undefined
        );
        bool b_is_backbone
        (
          atom_b.GetName() == "HA" ||
          atom_b.GetName() == "H" ||
          atom_b.GetName() == "CA" ||
          atom_b.GetName() == "CB" ||
          atom_b.GetName() == "N" ||
          atom_b.GetName() == "C" ||
          atom_b.GetName() == "HA2" ||
          atom_b.GetName() == "HA3" ||
          atom_b == biol::GetAtomTypes().e_Undefined
        );

        // change the atom type if it is not already a backbone atom
        if( !a_is_backbone)
        {
          atom_a = biol::GetAtomTypes().CB;
        }
        if( !b_is_backbone)
        {
          atom_b = biol::GetAtomTypes().CB;
        }

        // store the new data
        new_restraints.PushBack
        (
          util::ShPtr< restraint::AtomDistance>
          (
            new restraint::AtomDistance
            (
              util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
              (
                new restraint::LocatorCoordinatesHydrogen
                (
                  ( *itr)->GetData().First()->GetChainID(),
                  ( *itr)->GetData().First()->GetSeqID(),
                  atom_a.GetName()
                )
              ),
              util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
              (
                new restraint::LocatorCoordinatesHydrogen
                (
                  ( *itr)->GetData().Second()->GetChainID(),
                  ( *itr)->GetData().Second()->GetSeqID(),
                  atom_b.GetName()
                )
              ),
              distance
            )
          )
        );
      }
      return new_restraints;
    }

    //! @brief default constructor
    NmrFileConvert::NmrFileConvert() :
      m_InputFileFlag( new command::FlagStatic( "input_file", "file to which nmr data comes from")),
      m_InputFileParam
      (
        new command::Parameter
        (
          "input_path",
          "\tpath and name of the inputfile",
          ""
        )
      ),
      m_InputFileTypeParam
      (
        new command::Parameter
        (
          "input_type",
          "give the formatting style of the input file",
          command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "ROSETTA", "star")),
          "ROSETTA"
        )
      ),
      m_InputDataTypeParam
      (
        new command::Parameter
        (
          "nmr_type",
          "give the type of nmr data",
          command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "noe", "rdc")),
          "noe"
        )
      ),
      m_OutputFileFlag
      (
        new command::FlagStatic
        (
          "output_file",
          "file to which output should be sent"
        )
      ),
      m_OutputFileParam
      (
        new command::Parameter
        (
          "output_file",
          "\tpath and name of the output file",
          ""
        )
      ),
      m_OutputFileTypeParam
      (
        new command::Parameter
        (
          "output_type",
          "give the formatting style of the output file",
           command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "ROSETTA", "star")),
           "ROSETTA"
         )
       ),
      m_PdbFileFlag
      (
        new command::FlagStatic
        (
          "pdb_file",
          "file for determining the sequence offset",
          command::Parameter
          (
            "pdb_file",
            "\tpath and name of the pdb file",
            "[1AAJ.pdb]"
          )
        )
      ),
      m_RandomSubsetFlag
      (
        new command::FlagStatic
        (
          "random_subset",
          "get a random subset of the restraints given",
          command::Parameter
          (
            "nr_restraints",
            "\tnumber of restraints to get",
            util::Format()( util::GetUndefinedSize_t())
          )
        )
      ),
      m_SequenceSeparationFlag
      (
        new command::FlagStatic
        (
          "sequence_separation",
          "minimum sequence distance between two restraints",
          command::Parameter
          (
            "sequence_sep",
            "\tminimum sequence distance between two restraints",
            "0"
          )
        )
      ),
      m_GetCentroidDistancesFlag
      (
        new command::FlagStatic
        (
          "centroid_distances",
          "will convert side chain NOEs to centroid atom distances by adding 1 angstrom for every bond back to the backbone.  Only necessary in Rosetta files"
        )
      ),
      m_KBPotentialsPrefixFlag
      (
        new command::FlagStatic
        (
          "kb_potentials_prefix",
          "flag for passing KB potentials prefix, if given, ROSETTA output will use KB score rather than the default, bounded",
          command::Parameter
          (
            "kb_potentials_path",
            "location of KB potential files",
            ""
          )
        )
      ),
      m_RosettaWeightFlag
      (
        new command::FlagStatic
        (
          "rosetta_weight",
          "flag for setting weight for rosetta scores",
          command::Parameter
          (
            "rosetta_score_weight",
            "weight for rosetta scores",
            "1"
          )
        )
      ),
      m_AdjustSeqIDFlag
      (
        new command::FlagStatic( "adjust_seq_id", "adjust seq id based on starting PDB ID in protein model")
      )
      {
        m_InputFileFlag->PushBack( m_InputFileParam);
        m_InputFileFlag->PushBack( m_InputFileTypeParam);
        m_InputFileFlag->PushBack( m_InputDataTypeParam);
        m_OutputFileFlag->PushBack( m_OutputFileParam);
        m_OutputFileFlag->PushBack( m_OutputFileTypeParam);
      }

      const ApplicationType NmrFileConvert::NmrFileConvert_Instance
      (
        GetAppGroups().AddAppToGroup( new NmrFileConvert(), GetAppGroups().e_Restraint)
      );

  } // namespace app
} // namespace bcl
