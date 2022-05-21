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
#include "assemble/bcl_assemble_collector_common_aa.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CollectorCommonAA::s_Instance
    (
      util::Enumerated< find::CollectorInterface< storage::VectorND< 2, util::SiPtrList< const biol::AABase> >, storage::VectorND< 2, ProteinModel> > >::AddInstance
      (
        new CollectorCommonAA()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorCommonAA::CollectorCommonAA()
    {
    }

    //! @brief
    //! @param COLLECTOR ???
    CollectorCommonAA::CollectorCommonAA( const CollectorCommonAA &COLLECTOR)
    {
    }

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorSSEFurthest which is a copy of this
    CollectorCommonAA *CollectorCommonAA::Clone() const
    {
      return new CollectorCommonAA( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorCommonAA::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &CollectorCommonAA::GetAlias() const
    {
      static const std::string s_alias( "CollectorCommonAA");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorCommonAA::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Collects all amino acids that two proteins have in common.");

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief collect defined AAs in chains
    //! @param PROTEIN_MODEL model containing chains
    //! @return list of defined AAs
    util::SiPtrList< const biol::AABase> CollectorCommonAA::CollectDefinedAAsInChains
    (
      const ProteinModel &PROTEIN_MODEL
    )
    {
      // create list to hold the collected amino acids that have defined coordinates within the protein model
      util::SiPtrList< const biol::AABase> amino_acids;

      // loop over the chains of the two protein models
      for
      (
        storage::Vector< util::ShPtr< Chain> >::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // loop over the aa sequences of the chains denoted by "chain_itr_a" and "chain_itr_b"
        for
        (
          storage::Vector< util::ShPtr< biol::AABase> >::const_iterator
            aa_itr( ( *chain_itr)->GetSequence()->Begin()),
            aa_itr_end( ( *chain_itr)->GetSequence()->End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          // check if the coordinates are defined
          if
          (
            coord::AreDefinedCoordinates( ( *aa_itr)->GetAtomCoordinates())
          )
          {
            // add the amino acid denoted by "aa_itr" to first group in "amino_acids"
            amino_acids.PushBack( util::SiPtr< const biol::AABase>( **aa_itr));
          }
        }
      }

      // return the amino acids that are common to both proteins of "PROTEIN_MODELS"
      return amino_acids;
    }

    //! @brief collect defined SSEs in chains
    //! @param PROTEIN_MODEL model containing SSEs
    //! @return list of defined SSEs
    util::SiPtrList< const biol::AABase> CollectorCommonAA::CollectDefinedAAsInSSEs
    (
      const ProteinModel &PROTEIN_MODEL
    )
    {
      // create list to hold the collected amino acids that have defined coordinates within the protein model
      const util::SiPtrVector< const biol::AABase> amino_acids_initial( PROTEIN_MODEL.GetAminoAcids());
      util::SiPtrList< const biol::AABase> amino_acids;

      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          aa_itr( amino_acids_initial.Begin()),
          aa_itr_end( amino_acids_initial.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // check if the coordinates are defined
        if
        (
          coord::AreDefinedCoordinates( ( *aa_itr)->GetAtomCoordinates())
        )
        {
          // add the amino acid denoted by "aa_itr" to first group in "amino_acids"
          amino_acids.PushBack( *aa_itr);
        }
      }

      // return the amino acids that are common to both proteins of "PROTEIN_MODELS"
      return amino_acids;
    }

    //! @brief
    //! @param AAS_DEFINED ???
    //! @return
    storage::VectorND< 2, util::SiPtrList< const biol::AABase> >
    CollectorCommonAA::Collect
    (
      const storage::VectorND< 2, util::SiPtr< const util::SiPtrList< const biol::AABase> > > &AAS_DEFINED
    )
    {
      // create VectorND to hold the collected amino acids that are common to both protein models
      storage::VectorND< 2, util::SiPtrList< const biol::AABase> > amino_acids;

      util::SiPtrList< const biol::AABase>::const_iterator
        aa_itr_a( AAS_DEFINED.First()->Begin()),
        aa_itr_a_end( AAS_DEFINED.First()->End()),
        aa_itr_b( AAS_DEFINED.Second()->Begin()),
        aa_itr_b_end( AAS_DEFINED.Second()->End());

      while( aa_itr_a != aa_itr_a_end && aa_itr_b != aa_itr_b_end)
      {
        // match up chain ids
        while( aa_itr_a != aa_itr_a_end && ( *aa_itr_a)->GetChainID() < ( *aa_itr_b)->GetChainID())
        {
          ++aa_itr_a;
        }

        // break if there are no aas left in one of the lists
        if( aa_itr_a == aa_itr_a_end)
        {
          break;
        }

        while( aa_itr_b != aa_itr_b_end && ( *aa_itr_b)->GetChainID() < ( *aa_itr_a)->GetChainID())
        {
          ++aa_itr_b;
        }

        // break if there are no aas left in one of the lists
        if( aa_itr_b == aa_itr_b_end)
        {
          break;
        }

        // continue if they still do not match up
        if( ( *aa_itr_a)->GetChainID() != ( *aa_itr_b)->GetChainID())
        {
          continue;
        }

        // assert that the chain ids do indeed match up
        BCL_Assert
        (
            ( *aa_itr_a)->GetChainID() == ( *aa_itr_b)->GetChainID(),
            "chain id's do not match the first is :" + util::Format()( ( ( *aa_itr_a)->GetChainID()))
            + ": and the second is :"
            + util::Format()( ( ( *aa_itr_b)->GetChainID())) + ":"
        );

        // match up seq ids
        while( aa_itr_a != aa_itr_a_end && ( *aa_itr_a)->GetChainID() == ( *aa_itr_b)->GetChainID() && ( *aa_itr_a)->GetSeqID() < ( *aa_itr_b)->GetSeqID())
        {
          ++aa_itr_a;
        }

        // break if there are no aas left in one of the lists
        if( aa_itr_a == aa_itr_a_end)
        {
          break;
        }

        while( aa_itr_b != aa_itr_b_end && ( *aa_itr_a)->GetChainID() == ( *aa_itr_b)->GetChainID() && ( *aa_itr_b)->GetSeqID() < ( *aa_itr_a)->GetSeqID())
        {
          ++aa_itr_b;
        }

        // break if there are no aas left in one of the lists
        if( aa_itr_b == aa_itr_b_end)
        {
          break;
        }

        // continue if they still do not match up
        if( ( *aa_itr_a)->GetChainID() != ( *aa_itr_b)->GetChainID() || ( *aa_itr_a)->GetSeqID() != ( *aa_itr_b)->GetSeqID())
        {
          continue;
        }

        // assert that the seq ids do indeed match up
        BCL_Assert
        (
               ( *aa_itr_a)->GetChainID() == ( *aa_itr_b)->GetChainID()
            && ( *aa_itr_a)->GetSeqID()   == ( *aa_itr_b)->GetSeqID()
//            && ( *aa_itr_a)->GetPdbID()   == ( *aa_itr_b)->GetPdbID()
            ,
            "chain id's and seq ids should match... the first is chain :"
            + util::Format()( ( ( *aa_itr_a)->GetChainID())) + ": seq id :"
            + util::Format()( ( *aa_itr_a)->GetSeqID()) + ": pdb id :"
            + util::Format()( ( *aa_itr_a)->GetPdbID()) + ": \n"
            + "and the second is chainid :"
            + util::Format()( ( ( *aa_itr_b)->GetChainID())) + ": seq id :"
            + util::Format()( ( *aa_itr_b)->GetSeqID()) + ": pdb id :"
            + util::Format()( ( *aa_itr_b)->GetPdbID()) + ": \n"
        );

        // as long as amino acids do match, store them in amino_acids and got to next pair
        while( aa_itr_a != aa_itr_a_end && aa_itr_b != aa_itr_b_end && ( *aa_itr_a)->GetChainID() == ( *aa_itr_b)->GetChainID() && ( *aa_itr_a)->GetSeqID() == ( *aa_itr_b)->GetSeqID())
        {
          // make sure amino acids are of the same type
          BCL_Assert
          (
            ( **aa_itr_a) == ( **aa_itr_b),
            "amino acids have different type: "
            + util::Format()( ( *aa_itr_a)->GetType()) + " != " + util::Format()( ( *aa_itr_b)->GetType())
            + " first aa is PdbID" + util::Format()( ( *aa_itr_a)->GetPdbID()) + " in chain "
            + ( *aa_itr_a)->GetChainID()
            + " second aa is PdbID" + util::Format()( ( *aa_itr_b)->GetPdbID()) + " in chain "
            + ( *aa_itr_b)->GetChainID()
          );

          // store amino acids
          amino_acids.First().PushBack( *aa_itr_a);
          amino_acids.Second().PushBack( *aa_itr_b);
          ++aa_itr_a;
          ++aa_itr_b;
        }
      }

      return amino_acids;
    }

    //! @brief CollectCommonCoordinates is for aligning residues and giving corresponding aligned lists of atoms
    //! @param AAS_DEFINED the two list of amino acids which will be aligned to get their coordinates
    //! @param ATOM_TYPES the atom types in the residues for which you want the coordinates
    //! @return two lists of coordinates which are the aligned coordinates of the common aas and desired atoms
    storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> >
    CollectorCommonAA::CollectCommonCoordinates
    (
      const storage::VectorND< 2, util::SiPtr< const util::SiPtrList< const biol::AABase> > > &AAS_DEFINED,
      const storage::Set< biol::AtomType> &ATOM_TYPES,
      size_t &NR_COMMON_RESIDUES
    )
    {
      // create VectorND to hold the collected atom_coordinates that are common to both lists of AABases
      storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > all_atom_coordinates;

      size_t nr_common_residues( 0);

      util::SiPtrList< const biol::AABase>::const_iterator
        aa_itr_a( AAS_DEFINED.First()->Begin()),
        aa_itr_a_end( AAS_DEFINED.First()->End()),
        aa_itr_b( AAS_DEFINED.Second()->Begin()),
        aa_itr_b_end( AAS_DEFINED.Second()->End());

      while( aa_itr_a != aa_itr_a_end && aa_itr_b != aa_itr_b_end)
      {
        // try to match up chain ids by  moving "aa_itr_a"
        while( aa_itr_a != aa_itr_a_end && ( *aa_itr_a)->GetChainID() < ( *aa_itr_b)->GetChainID())
        {
          ++aa_itr_a;
        }

        // break if there are no aas left in list a
        if( aa_itr_a == aa_itr_a_end)
        {
          break;
        }

        // try to match up chain ids by  moving "aa_itr_b"
        while( aa_itr_b != aa_itr_b_end && ( *aa_itr_b)->GetChainID() < ( *aa_itr_a)->GetChainID())
        {
          ++aa_itr_b;
        }

        // break if there are no aas left in list b
        if( aa_itr_b == aa_itr_b_end)
        {
          break;
        }

        // continue if they still do not match up i.e. one whole chain is missing in one of the lists of aas
        if( ( *aa_itr_a)->GetChainID() != ( *aa_itr_b)->GetChainID())
        {
          continue;
        }

        // assert that the chain ids do indeed match up
        BCL_Assert
        (
            ( *aa_itr_a)->GetChainID() == ( *aa_itr_b)->GetChainID(),
            "chain id's do not match the first is :" + util::Format()( ( ( *aa_itr_a)->GetChainID()))
            + ": and the second is :"
            + util::Format()( ( ( *aa_itr_b)->GetChainID())) + ":"
        );

        // match up seq ids by trying to move "aa_itr_a"
        while
        (
          aa_itr_a != aa_itr_a_end &&
          ( *aa_itr_a)->GetChainID() == ( *aa_itr_b)->GetChainID() &&
          ( *aa_itr_a)->GetSeqID() < ( *aa_itr_b)->GetSeqID()
        )
        {
          ++aa_itr_a;
        }

        // break if there are no aas left in list a
        if( aa_itr_a == aa_itr_a_end)
        {
          break;
        }

        // match up seq ids by trying to move "aa_itr_b"
        while
        (
          aa_itr_b != aa_itr_b_end &&
          ( *aa_itr_a)->GetChainID() == ( *aa_itr_b)->GetChainID() &&
          ( *aa_itr_b)->GetSeqID() < ( *aa_itr_a)->GetSeqID()
        )
        {
          ++aa_itr_b;
        }

        // break if there are no aas left in list b
        if( aa_itr_b == aa_itr_b_end)
        {
          break;
        }

        // continue if they still do not match up
        if
        (
          ( *aa_itr_a)->GetChainID() != ( *aa_itr_b)->GetChainID() ||
          ( *aa_itr_a)->GetSeqID() != ( *aa_itr_b)->GetSeqID()
        )
        {
          continue;
        }

        // assert that the seq ids do indeed match up
        BCL_Assert
        (
               ( *aa_itr_a)->GetChainID() == ( *aa_itr_b)->GetChainID()
            && ( *aa_itr_a)->GetSeqID()   == ( *aa_itr_b)->GetSeqID()
//            && ( *aa_itr_a)->GetPdbID()   == ( *aa_itr_b)->GetPdbID()
            ,
            "chain id's and seq ids should match... the first is chain :"
            + util::Format()( ( ( *aa_itr_a)->GetChainID())) + ": seq id :"
            + util::Format()( ( *aa_itr_a)->GetSeqID()) + ": pdb id :"
            + util::Format()( ( *aa_itr_a)->GetPdbID()) + ": \n"
            + "and the second is chainid :"
            + util::Format()( ( ( *aa_itr_b)->GetChainID())) + ": seq id :"
            + util::Format()( ( *aa_itr_b)->GetSeqID()) + ": pdb id :"
            + util::Format()( ( *aa_itr_b)->GetPdbID()) + ": \n"
        );

        // as long as amino acids do match, store them in amino_acids and got to next pair
        while
        (
          aa_itr_a != aa_itr_a_end &&
          aa_itr_b != aa_itr_b_end &&
          ( *aa_itr_a)->GetChainID() == ( *aa_itr_b)->GetChainID() &&
          ( *aa_itr_a)->GetSeqID() == ( *aa_itr_b)->GetSeqID()
        )
        {
          // make sure amino acids are of the same type
          BCL_Assert
          (
            ( **aa_itr_a) == ( **aa_itr_b),
            "amino acids have different type: "
            + util::Format()( ( *aa_itr_a)->GetType()) + " != " + util::Format()( ( *aa_itr_b)->GetType())
            + " first aa is PdbID" + util::Format()( ( *aa_itr_a)->GetPdbID()) + " in chain "
            + ( *aa_itr_a)->GetChainID()
            + " second aa is PdbID" + util::Format()( ( *aa_itr_b)->GetPdbID()) + " in chain "
            + ( *aa_itr_b)->GetChainID()
          );

          // store coordinates
          const util::SiPtrVector< const linal::Vector3D> temp_coords_a( ( *aa_itr_a)->GetAtomCoordinates( ATOM_TYPES));
          const util::SiPtrVector< const linal::Vector3D> temp_coords_b( ( *aa_itr_b)->GetAtomCoordinates( ATOM_TYPES));

          // skip undefined coordinates
          if( !coord::AreDefinedCoordinates( temp_coords_a))
          {
            BCL_MessageVrb
            (
              "skipping aa with undefined coordinates: chain" + util::Format()( ( *aa_itr_a)->GetChainID()) + " pdbid: " +
              util::Format()( ( *aa_itr_a)->GetPdbID())
            );
            ++aa_itr_a;
            ++aa_itr_b;
            continue;
          }
          if( !coord::AreDefinedCoordinates( temp_coords_b))
          {
            BCL_MessageVrb
            (
              "skipping aa with undefined coordinates: chain" + util::Format()( ( *aa_itr_b)->GetChainID()) + " pdbid: " +
              util::Format()( ( *aa_itr_b)->GetPdbID())
            );
            ++aa_itr_a;
            ++aa_itr_b;
            continue;
          }

          // add coordinates to the end of appropriate SiPtrVectors of "COMMON_AAS"
          all_atom_coordinates.First().Append( temp_coords_a);
          all_atom_coordinates.Second().Append( temp_coords_b);
          ++nr_common_residues;
          ++aa_itr_a;
          ++aa_itr_b;
        }
      }

      NR_COMMON_RESIDUES = nr_common_residues;
      return all_atom_coordinates;
    }

    //! @brief CollectCommonCoordinates is for aligning residues and giving corresponding aligned lists of atoms
    //! @param AAS_DEFINED_A the list of amino acids which will be aligned to get their coordinates from first model
    //! @param AAS_DEFINED_B the list of amino acids which will be aligned to get their coordinates from second model
    //! @param ATOM_TYPES the atom types in the residues for which you want the coordinates
    //! @param NR_COMMON_RESIDUES number of common residues, this value will be overwritten by the function
    //! @return two lists of coordinates which are the aligned coordinates of the common aas and desired atoms
    storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> >
    CollectorCommonAA::CollectCommonCoordinates
    (
      const util::SiPtrList< const biol::AABase> &AAS_DEFINED_A,
      const util::SiPtrList< const biol::AABase> &AAS_DEFINED_B,
      const storage::Set< biol::AtomType> &ATOM_TYPES,
      size_t &NR_COMMON_RESIDUES
    )
    {
      return CollectCommonCoordinates
      (
        storage::VectorND< 2, util::SiPtr< const util::SiPtrList< const biol::AABase> > >
        (
          util::ToSiPtr( AAS_DEFINED_A),
          util::ToSiPtr( AAS_DEFINED_B)
        ),
        ATOM_TYPES,
        NR_COMMON_RESIDUES
      );

    }

    //! @brief ???
    //! @param PROTEIN_MODELS ???
    //! @return returns
    storage::VectorND< 2, util::SiPtrList< const biol::AABase> >
    CollectorCommonAA::Collect( const storage::VectorND< 2, ProteinModel> &PROTEIN_MODELS) const
    {
      // create VectorND "" to hold the collected amino acids that are common to both protein models
      storage::VectorND< 2, util::SiPtrList< const biol::AABase> > amino_acids;

      // loop over the chains of the two protein models
      for
      (
        storage::Vector< util::ShPtr< Chain> >::const_iterator
          chain_itr_a( PROTEIN_MODELS.First().GetChains().Begin()),
          chain_itr_a_end( PROTEIN_MODELS.First().GetChains().End()),
          chain_itr_b( PROTEIN_MODELS.Second().GetChains().Begin()),
          chain_itr_b_end( PROTEIN_MODELS.Second().GetChains().End());
        chain_itr_a != chain_itr_a_end && chain_itr_b != chain_itr_b_end;
        ++chain_itr_a, ++chain_itr_b
      )
      {
        // make sure chain id's match
        BCL_Assert( ( *chain_itr_a)->GetChainID() == ( *chain_itr_b)->GetChainID(), "chain ids do not match");

        // loop over the aa sequences of the chains denoted by "chain_itr_a" and "chain_itr_b"
        for
        (
          storage::Vector< util::ShPtr< biol::AABase> >::const_iterator
            aa_itr_a( ( *chain_itr_a)->GetSequence()->Begin()),
            aa_itr_a_end( ( *chain_itr_a)->GetSequence()->End()),
            aa_itr_b( ( *chain_itr_b)->GetSequence()->Begin()),
            aa_itr_b_end( ( *chain_itr_b)->GetSequence()->End());
          aa_itr_a != aa_itr_a_end && aa_itr_b != aa_itr_b_end;
          ++aa_itr_a, ++aa_itr_b
        )
        {
          // make sure amino acids are of the same type
          BCL_Assert
          (
            ( **aa_itr_a) == ( **aa_itr_b),
            "amino acids have different type: "
            + util::Format()( ( *aa_itr_a)->GetType()) + " != " + util::Format()( ( *aa_itr_b)->GetType())
          );

          // check if the coordinates are defined
          if
          (
               coord::AreDefinedCoordinates( ( *aa_itr_a)->GetAtomCoordinates())
            && coord::AreDefinedCoordinates( ( *aa_itr_b)->GetAtomCoordinates())
          )
          {
            // add the amino acid denoted by "aa_itr_a" to first group in "amino_acids"
            amino_acids.First().PushBack( util::ToSiPtr( **aa_itr_a));

            // add the amino acid denoted by "aa_itr_b" to first group in "amino_acids"
            amino_acids.Second().PushBack( util::ToSiPtr( **aa_itr_b));
          }
        }
      }

      // return the amino acids that are common to both proteins of "PROTEIN_MODELS"
      return amino_acids;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace assemble
} // namespace bcl
