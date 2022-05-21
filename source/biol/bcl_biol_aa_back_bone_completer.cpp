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
#include "biol/bcl_biol_aa_back_bone_completer.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_classes.h"
#include "biol/bcl_biol_aa_complete.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! the bond angle between the CA->C->N atoms where N is of the residue following the CA and C atom residue
    //! the value is 117.2 degrees (2.0455259 radians) and comes from reference
    //! Berkholz, D. et al. "Conformation Dependence of Backbone Geometry in Proteins" Structure 17( 2009): 1316-1325
    const double AABackBoneCompleter::s_CACNBondAngle( 2.0455259);

    //! the bond angle between the O->C->N atoms where N  is of the residue following the CA and C atom residue
    const double AABackBoneCompleter::s_OCNAngle( math::g_Pi - 0.5 * s_CACNBondAngle);

    //! the bond angle between the C->N->CA atoms where C is of the residue preceding the N and CA atom residue
    //! the value is 121.7 degrees (2.1240657 radians) and comes from reference
    //! Berkholz, D. et al. "Conformation Dependence of Backbone Geometry in Proteins" Structure 17( 2009): 1316-1325
    const double AABackBoneCompleter::s_CNCABondAngle( 2.1240657);

    //! the bond angle between the C->N->H atoms where C is of the residue preceding the N and H atom residue
    const double AABackBoneCompleter::s_CNHBondAngle( math::g_Pi - 0.5 * s_CNCABondAngle);

    //! the ideal omega bond angle CA->C->N->CA
    const double AABackBoneCompleter::s_OmegaCACNCA( math::g_Pi);

    //! the ideal omega bond angle CA->C->N->H
    const double AABackBoneCompleter::s_OmegaCACNH( 0.0);

    //! the ideal omega bond angle O->C->N->CA
    const double AABackBoneCompleter::s_OmegaOCNCA( 0.0);

    //! the ideal omega bond angle O->C->N->H
    const double AABackBoneCompleter::s_OmegaOCNH( math::g_Pi);

    // instantiate static instance
    const util::SiPtr< const util::ObjectInterface> AABackBoneCompleter::s_Instance
    (
      GetObjectInstances().AddInstance( new AABackBoneCompleter())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from two booleans whether or not to add missing amide hydrogen and carbonyl oxygens
    //! @param ADD_AMIDE_HYDROGENS whether or not to add missing amide hydrogens
    //! @param ADD_CARBONYL_OXYGENS whether or not to add missing carbonyl oxygens
    //! @param ADD_HA_HYDROGENS whether or not to add missing HA hydrogens
    AABackBoneCompleter::AABackBoneCompleter
    (
      const bool ADD_AMIDE_HYDROGENS,
      const bool ADD_CARBONYL_OXYGENS,
      const bool ADD_HA_HYDROGENS
    ) :
      m_AddAmideHydrogens( ADD_AMIDE_HYDROGENS),
      m_AddCarbonylOxygens( ADD_CARBONYL_OXYGENS),
      m_AddHAHydrogens( ADD_HA_HYDROGENS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AABackBoneCompleter
    AABackBoneCompleter *AABackBoneCompleter::Clone() const
    {
      return new AABackBoneCompleter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AABackBoneCompleter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to complete backbones for a given ProteinModel
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return ProteinModel with backbone completed and of type AAComplete
    //! @return ShPtr to a new ProteinModel with completed backbone
    util::ShPtr< assemble::ProteinModel> AABackBoneCompleter::CompleteProteinModel
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // construct store for new chains to be constructed
      util::ShPtrVector< assemble::Chain> new_chains;

      // iterate over chains in PROTEIN_MODEL
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_end_itr( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_end_itr; ++chain_itr
      )
      {
        // pushback copy of this chain
        new_chains.PushBack( CompleteChain( **chain_itr));
      }

      // construct and return a new ProteinModel
      util::ShPtr< assemble::ProteinModel> new_model( new assemble::ProteinModel( new_chains));
      new_model->SetProteinModelData( PROTEIN_MODEL.GetProteinModelData());
      return new_model;
    }

    //! @brief function to complete backbones for a given Chain
    //! @param CHAIN Chain of interest
    //! @return Chain with backbone completed and of type AAComplete
    //! @return ShPtr to a new Chain with completed backbone
    util::ShPtr< assemble::Chain> AABackBoneCompleter::CompleteChain
    (
      const assemble::Chain &CHAIN
    ) const
    {
      // create storage for new SSEs
      util::ShPtrVector< assemble::SSE> new_sses;

      // iterate over SSEs in this chain
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( CHAIN.GetData().Begin()), sse_itr_end( CHAIN.GetData().End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        util::SiPtr< const AABase> sp_prev_aa;
        if( sse_itr != CHAIN.GetData().Begin())
        {
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr_prev( sse_itr);
          --sse_itr_prev;

          if( ( *sse_itr_prev)->GetLastAA()->DoesPrecede( *( *sse_itr)->GetFirstAA()))
          {
            sp_prev_aa = ( *sse_itr_prev)->GetLastAA();
          }
        }
        util::SiPtr< const AABase> sp_next_aa;
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr_next( sse_itr);
        ++sse_itr_next;
        if( sse_itr_next != sse_itr_end)
        {
          if( ( *sse_itr)->GetLastAA()->DoesPrecede( *( *sse_itr_next)->GetFirstAA()))
          {
            sp_next_aa = ( *sse_itr_next)->GetFirstAA();
          }
        }

        // construct a new SSE with complete backbone
        util::ShPtr< assemble::SSE> new_sse
        (
          new assemble::SSE( *CompleteAASequence( **sse_itr, sp_prev_aa, sp_next_aa), ( *sse_itr)->GetType())
        );
        // insert it into new_sses
        new_sses.PushBack( new_sse);
      }

      // construct a new chain with the sequence and the new sses and return it
      return util::ShPtr< assemble::Chain>( new assemble::Chain( CHAIN.GetSequence(), new_sses));
    }

    //! @brief function to complete backbones for a given AASequence
    //! @param AA_SEQUENCE AASequence of interest
    //! @param SP_PREV_AA amino acid previous to the sequence
    //! @param SP_NEXT_AA amino acid following the sequence
    //! @return AASequence with backbone completed and of type AAComplete
    //! @return ShPtr to a new AASequence with completed backbone
    util::ShPtr< AASequence> AABackBoneCompleter::CompleteAASequence
    (
      const AASequence &AA_SEQUENCE,
      const util::SiPtr< const AABase> &SP_PREV_AA,
      const util::SiPtr< const AABase> &SP_NEXT_AA
    ) const
    {
      // create storage for new amino acids
      util::ShPtrVector< AABase> new_amino_acids;

      // iterate over amino acids
      for
      (
        AASequence::const_iterator
          aa_itr( AA_SEQUENCE.Begin()), aa_itr_end( AA_SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        util::SiPtr< const AABase> sp_prev_aa;
        util::SiPtr< const AABase> sp_next_aa;
        if( aa_itr == AA_SEQUENCE.Begin())
        {
          sp_prev_aa = SP_PREV_AA;
        }
        else
        {
          sp_prev_aa = *( aa_itr - 1);
        }

        if( aa_itr + 1 == aa_itr_end)
        {
          sp_next_aa = SP_NEXT_AA;
        }
        else
        {
          sp_next_aa = *( aa_itr + 1);
        }

        // create a reference on the amino acid
        const AABase &amino_acid( **aa_itr);

        BCL_MessageDbg( "amino acid: " + amino_acid.GetIdentification());
        // make sure that this amino acid is back bone or complete
        BCL_Assert
        (
          amino_acid.GetAAClass() == GetAAClasses().e_AABackBone ||
          amino_acid.GetAAClass() == GetAAClasses().e_AAComplete,
          "To complete backbone the amino acid has to be either e_AABackBone or e_AAComplete not " +
          amino_acid.GetAAClass().GetName()
        );

        // create a new AAComplete with all the atoms
        util::ShPtr< AABase> new_amino_acid( new AAComplete( amino_acid));

        // search for the oxygen
        Atom oxygen( amino_acid.GetAtom( GetAtomTypes().O));

        // if oxygen is to be added and
        // it is not the last amino acid and
        // there is no oxygen already in the amino acid or it has undefined coordinates
        if
        (
          m_AddCarbonylOxygens &&
          sp_next_aa.IsDefined() &&
          !( oxygen.GetType().IsDefined() && oxygen.GetCoordinates().IsDefined())
        )
        {
          // generate a oxygen atom using this amino acid and N coordinates from next amino acid
          oxygen =
            GenerateOxygen
            (
              amino_acid,
              sp_next_aa->GetAtom( GetAtomTypes().N).GetCoordinates()
            );

          BCL_MessageDbg( "adding oxygen with coord " + util::Format()( oxygen.GetCoordinates()));

          // if the oxygen was correctly constructed
          if( oxygen.GetCoordinates().IsDefined())
          {
            // add it to new amino acid
            new_amino_acid->SetAtom( oxygen);
          }
        }

        // skip prolines, which don't have H
        if( amino_acid.GetType() != GetAATypes().PRO)
        {
          // search for the hydrogen
          Atom hydrogen( amino_acid.GetAtom( GetAtomTypes().H));

          // if hydrogen is to be added and
          // it is not first amino acid and
          // there is no hydrogen already in the amino acid or it has undefined coordinates
          if
          (
            m_AddAmideHydrogens &&
            sp_prev_aa.IsDefined() &&
            !( hydrogen.GetType().IsDefined() && hydrogen.GetCoordinates().IsDefined())
          )
          {
            // generate a hydrogen atom using this amino acid and C coordinates from previous amino acid
            hydrogen =
              GenerateHydrogen
              (
                amino_acid,
                sp_prev_aa->GetAtom( GetAtomTypes().C).GetCoordinates()
              );

            BCL_MessageDbg
            (
              "adding amide hydrogen with coord " + util::Format()( hydrogen.GetCoordinates())
            );

            // if hydrogen was correctly constructed
            if( hydrogen.GetCoordinates().IsDefined())
            {
              // add it to new amino acid
              new_amino_acid->SetAtom( hydrogen);
            }
          }
        }

        // search for the HA (or HA3 for GLY)
        const bool is_glycine( amino_acid.GetType() == GetAATypes().GLY);
        Atom ha
        (
          is_glycine ?
            amino_acid.GetAtom( GetAtomTypes().HA3) :
            amino_acid.GetAtom( GetAtomTypes().HA)
        );

        // if HA is to be added and there is no hydrogen already in the amino acid or it has undefined coordinates
        if( m_AddHAHydrogens && !( ha.GetType().IsDefined() && ha.GetCoordinates().IsDefined()))
        {
          ha = GenerateHA( amino_acid);

          BCL_MessageDbg
          (
            "adding alpha hydrogen with coord " + util::Format()( ha.GetCoordinates())
          );

          bool determined_ha( ha.GetCoordinates().IsDefined());
          // if ha was correctly constructed
          if( determined_ha)
          {
            // add it to new amino acid
            new_amino_acid->SetAtom( ha);
          }

          // if this is a glycine
          if( is_glycine)
          {
            // get the HA2
            const Atom &ha2( amino_acid.GetAtom( GetAtomTypes().HA2));

            // if HA2 is not already present with defined coordinates
            if( !ha2.GetType().IsDefined() || !ha2.GetCoordinates().IsDefined())
            {
              if( determined_ha)
              {
                // calculate the coordinates of where HA2 should go
                const linal::Vector3D ha2_coordinates
                (
                  linal::CoordinatesTetrahedral
                  (
                    new_amino_acid->GetAtom( GetAtomTypes().CA).GetCoordinates(),
                    new_amino_acid->GetAtom( GetAtomTypes().N).GetCoordinates(),
                    new_amino_acid->GetAtom( GetAtomTypes().C).GetCoordinates(),
                    new_amino_acid->GetAtom( GetAtomTypes().HA3).GetCoordinates(),
                    GetAtomTypes().CA->GetBondLength( GetAtomTypes().HA)
                  )
                );

                // if the coordinates were found for HA2
                if( ha2_coordinates.IsDefined())
                {
                  // add it to new amino acid
                  new_amino_acid->SetAtom( Atom( ha2_coordinates, GetAtomTypes().HA2));
                }
              }
              else
              {
                const double bond_length( GetAtomTypes().CA->GetBondLength( GetAtomTypes().HA));
                const linal::Vector3D &position( new_amino_acid->GetAtom( GetAtomTypes().CA).GetCoordinates());
                // helper coordinates
                linal::Vector3D foot_point
                (
                  linal::CoordinatesTrigonal
                  (
                    position,
                    new_amino_acid->GetAtom( GetAtomTypes().N).GetCoordinates(),
                    new_amino_acid->GetAtom( GetAtomTypes().C).GetCoordinates(),
                    bond_length * std::cos( 54.75 / 180 * math::g_Pi)
                  )
                );

                linal::Vector3D offset
                (
                  bond_length * std::sin( 54.75 / 180 * math::g_Pi) *
                  linal::CrossProduct
                  (
                    new_amino_acid->GetAtom( GetAtomTypes().N).GetCoordinates() - position,
                    new_amino_acid->GetAtom( GetAtomTypes().C).GetCoordinates() - position
                  ).Normalize()
                );

                new_amino_acid->SetAtom( Atom( foot_point + offset, GetAtomTypes().HA));
                new_amino_acid->SetAtom( Atom( foot_point - offset, GetAtomTypes().HA2));
              }
            }
          } // is glycine
        } // add HA

        // create a new amino acid from this amino acid with side chains
        new_amino_acids.PushBack( new_amino_acid);
      }

      // construct a new Sequence and return it
      return util::ShPtr< AASequence>
      (
        new AASequence( new_amino_acids, AA_SEQUENCE.GetChainID(), AA_SEQUENCE.GetFastaHeader())
      );
    }

    //! @brief function to generate a missing hydrogen from amino acid and previous C coordinates
    //! @param AMINO_ACID amino acid of interest
    //! @param PREV_C_COORD Coordinates of C of previous amino acid
    //! @param BOND_LENGTH N-H bond length
    //! @return generated hydrogen atom
    Atom AABackBoneCompleter::GenerateHydrogen
    (
      const AABase &AMINO_ACID,
      const linal::Vector3D &PREV_C_COORD,
      const double BOND_LENGTH
    )
    {
      // calculate the coordinates
      const linal::Vector3D hdyrogen_coordinates
      (
        linal::CoordinatesTrigonal
        (
          AMINO_ACID.GetAtom( GetAtomTypes().N).GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().CA).GetCoordinates(),
          PREV_C_COORD,
          BOND_LENGTH
        )
      );

      // construct and return hydrogen atom
      return Atom( hdyrogen_coordinates, GetAtomTypes().H);
    }

    //! @brief function to generate a missing oxygen from a amino acid and next N coordinates
    //! @param AMINO_ACID amino acid of interest
    //! @param NEXT_N_COORD Coordinates of N of next amino acid
    //! @param BOND_LENGTH C=O bond length
    //! @return generated oxygen atom
    Atom AABackBoneCompleter::GenerateOxygen
    (
      const AABase &AMINO_ACID,
      const linal::Vector3D &NEXT_N_COORD,
      const double BOND_LENGTH
    )
    {
      // calculate the coordinates
      const linal::Vector3D oxygen_coordinates
      (
        linal::CoordinatesTrigonal
        (
          AMINO_ACID.GetAtom( GetAtomTypes().C).GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().CA).GetCoordinates(),
          NEXT_N_COORD,
          BOND_LENGTH
        )
      );

      // construct and return oxygen atom
      return Atom( oxygen_coordinates, GetAtomTypes().O);
    }

    //! @brief function to generate a missing HA from a amino acid
    //! @param AMINO_ACID amino acid of interest
    //! @param BOND_LENGTH CA-HA bond length
    //! @return generated oxygen atom
    Atom AABackBoneCompleter::GenerateHA
    (
      const AABase &AMINO_ACID,
      const double BOND_LENGTH
    )
    {
      // calculate the coordinates
      const linal::Vector3D ha_coordinates
      (
        linal::CoordinatesTetrahedral
        (
          AMINO_ACID.GetAtom( GetAtomTypes().CA).GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().N).GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().C).GetCoordinates(),
          AMINO_ACID.GetFirstSidechainAtom().GetCoordinates(),
          BOND_LENGTH
        )
      );

      // construct and return HA (or HA3) atom
      return Atom
      (
        ha_coordinates,
        AMINO_ACID.GetType() == GetAATypes().GLY ? GetAtomTypes().HA3 : GetAtomTypes().HA
      );
    }

    //! @brief function to generate the N of the neighboring amino acid in the peptide bond
    //! @param AMINO_ACID amino acid of interest, with psi angle defined (using oxgyen)
    //! @param BOND_LENGTH C-N bond length
    //! @return generated nitrogen (c-terminal direction)
    Atom AABackBoneCompleter::GenerateN
    (
      const AABase &AMINO_ACID,
      const double BOND_LENGTH
    )
    {
      // calculate the coordinates
      const linal::Vector3D n_coordinates
      (
        linal::CoordinatesDihedral
        (
          AMINO_ACID.GetAtom( GetAtomTypes().C).GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().CA).GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().N).GetCoordinates(),
          BOND_LENGTH,
          s_CACNBondAngle,
          AMINO_ACID.Psi()
        )
      );

      // construct and return N atom
      return Atom( n_coordinates, GetAtomTypes().N);
    }

    //! @brief function to generate the C of the neighboring amino acid in the petide bond
    //! @param AMINO_ACID amino acid of interest
    //! @param PHI the phi angle
    //! @param BOND_LENGTH the C-N bond length
    //! @return generated carbon (n-terminal direction)
    Atom AABackBoneCompleter::GenerateC
    (
      const AABase &AMINO_ACID,
      const double PHI,
      const double BOND_LENGTH
    )
    {
      const linal::Vector3D desired_c_coordinates
      (
        linal::CoordinatesDihedral
        (
          AMINO_ACID.GetAtom( GetAtomTypes().N).GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().CA).GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().C).GetCoordinates(),
          BOND_LENGTH,
          s_CNCABondAngle,
          PHI
        )
      );

      // end
      return Atom( desired_c_coordinates, GetAtomTypes().C);
    }

    //! @brief function to generate C, O and CA on an N-terminal amino acid
    //! @param AMINO_ACID amino acid of interest
    //! @param PHI the phi angle on the AMINO acid
    //! @return map of atoms with C, O and CA atoms, peptide bond to AMINO_ACID
    storage::Map< AtomType, Atom> AABackBoneCompleter::GenerateCOCA
    (
      const AABase &AMINO_ACID,
      const double PHI
    )
    {
      const Atom carbon( GenerateC( AMINO_ACID, PHI));
      const Atom carbon_alpha
      (
        linal::CoordinatesDihedral
        (
          carbon.GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().N).GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().CA).GetCoordinates(),
          GetAtomTypes().C->GetBondLength( GetAtomTypes().CA),
          s_CACNBondAngle, // trigonal
          s_OmegaCACNCA    // omega
        ),
        GetAtomTypes().CA
      );
      const Atom oxygen
      (
        linal::CoordinatesTrigonal
        (
          carbon.GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().N).GetCoordinates(),
          carbon_alpha.GetCoordinates(),
          GetAtomTypes().C->GetBondLength( GetAtomTypes().O)
        ),
        GetAtomTypes().O
      );

      storage::Map< AtomType, Atom> atoms;
      atoms[ carbon.GetType()] = carbon;
      atoms[ oxygen.GetType()] = oxygen;
      atoms[ carbon_alpha.GetType()] = carbon_alpha;

      // end
      return atoms;
    }

    //! @brief function to generate H, N and CA on an C-terminal amino acid
    //! @param AMINO_ACID amino acid of interest
    //! @return map of atoms with H, N and CA atoms, peptide bond to AMINO_ACID
    storage::Map< AtomType, Atom> AABackBoneCompleter::GenerateHNCA
    (
      const AABase &AMINO_ACID
    )
    {
      const Atom nitrogen( GenerateN( AMINO_ACID));
      const Atom carbon_alpha
      (
        linal::CoordinatesDihedral
        (
          nitrogen.GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().C).GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().O).GetCoordinates(),
          GetAtomTypes().N->GetBondLength( GetAtomTypes().CA), // bond length
          s_CNCABondAngle, // trigonal
          s_OmegaOCNCA     // omega
        ),
        GetAtomTypes().CA
      );
      const Atom hydrogen
      (
        linal::CoordinatesTrigonal
        (
          nitrogen.GetCoordinates(),
          AMINO_ACID.GetAtom( GetAtomTypes().C).GetCoordinates(),
          carbon_alpha.GetCoordinates(),
          GetAtomTypes().N->GetBondLength( GetAtomTypes().H)
        ),
        GetAtomTypes().H
      );

      // all atoms
      storage::Map< AtomType, Atom> atoms;
      atoms[ hydrogen.GetType()] = hydrogen;
      atoms[ nitrogen.GetType()] = nitrogen;
      atoms[ carbon_alpha.GetType()] = carbon_alpha;

      // end
      return atoms;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AABackBoneCompleter::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AddAmideHydrogens, ISTREAM);
      io::Serialize::Read( m_AddCarbonylOxygens, ISTREAM);
      io::Serialize::Read( m_AddHAHydrogens, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AABackBoneCompleter::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AddAmideHydrogens, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AddCarbonylOxygens, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AddHAHydrogens, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
