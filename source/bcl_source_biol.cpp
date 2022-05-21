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
#include "biol/bcl_biol_aa_back_bone.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_classes.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AABackBone::s_Instance
    (
      GetObjectInstances().AddInstance( new AABackBone())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AABackBone::AABackBone() :
      AABase(),
      m_N( GetAtomTypes().N),
      m_CA( GetAtomTypes().CA),
      m_C( GetAtomTypes().C),
      m_O( GetAtomTypes().O),
      m_FirstSidechainAtom(),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_N, m_CA, m_C, m_O, m_FirstSidechainAtom))
    {
    }

    //! @brief construct AABackBone from util::ShPtr to AAData
    //! @brief SP_AA_DATA ShPTr to AAData to be copied
    AABackBone::AABackBone( const util::ShPtr< AAData> &SP_AA_DATA) :
      AABase( SP_AA_DATA),
      m_N( GetAtomTypes().N),
      m_CA( GetAtomTypes().CA),
      m_C( GetAtomTypes().C),
      m_O( GetAtomTypes().O),
      m_FirstSidechainAtom( GetType()->GetFirstSidechainAtomType()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_N, m_CA, m_C, m_O, m_FirstSidechainAtom))
    {
    }

    //! @brief constructor from AABase
    //! @param AA_BASE AABase
    AABackBone::AABackBone( const AABase &AA_BASE) :
      AABase( AA_BASE),
      m_N( GetAtomTypes().N),
      m_CA( GetAtomTypes().CA),
      m_C( GetAtomTypes().C),
      m_O( GetAtomTypes().O),
      m_FirstSidechainAtom( GetType()->GetFirstSidechainAtomType()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_N, m_CA, m_C, m_O, m_FirstSidechainAtom))
    {
    }

    //! @brief construct AACaCb from AABASE and util::SiPtrVector< const Atom> ATOMS
    //! @param AA_BASE AABase
    //! @param ATOMS SiPtrVector of Atoms
    AABackBone::AABackBone( const AABase &AA_BASE, const util::SiPtrVector< const Atom> &ATOMS) :
      AABase( AA_BASE),
      m_N(  Atom::FindAtom( ATOMS, GetAtomTypes().N)->GetType().IsDefined()  ? *Atom::FindAtom( ATOMS, GetAtomTypes().N)  : GetAtomTypes().N),
      m_CA( Atom::FindAtom( ATOMS, GetAtomTypes().CA)->GetType().IsDefined() ? *Atom::FindAtom( ATOMS, GetAtomTypes().CA) : GetAtomTypes().CA),
      m_C(  Atom::FindAtom( ATOMS, GetAtomTypes().C)->GetType().IsDefined()  ? *Atom::FindAtom( ATOMS, GetAtomTypes().C)  : GetAtomTypes().C),
      m_O(  Atom::FindAtom( ATOMS, GetAtomTypes().O)->GetType().IsDefined()  ? *Atom::FindAtom( ATOMS, GetAtomTypes().O)  : GetAtomTypes().O),
      m_FirstSidechainAtom( Atom::FindAtom( ATOMS, GetType()->GetFirstSidechainAtomType())->GetType().IsDefined() ? *Atom::FindAtom( ATOMS, GetType()->GetFirstSidechainAtomType()) : GetType()->GetFirstSidechainAtomType()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_N, m_CA, m_C, m_O, m_FirstSidechainAtom))
    {
    }

    //! @brief copy constructor
    //! @param AA_BACKBONE AABackBone to be copied
    AABackBone::AABackBone( const AABackBone &AA_BACKBONE) :
      AABase( AA_BACKBONE),
      m_N( AA_BACKBONE.m_N),
      m_CA( AA_BACKBONE.m_CA),
      m_C( AA_BACKBONE.m_C),
      m_O( AA_BACKBONE.m_O),
      m_FirstSidechainAtom( AA_BACKBONE.m_FirstSidechainAtom),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_N, m_CA, m_C, m_O, m_FirstSidechainAtom))
    {
    }

    //! @brief constructor from CaCb amino acids
    //! @param AA_CACB AACaCb to be copied
    AABackBone::AABackBone( const AACaCb &AA_CACB) :
      AABase( AA_CACB),
      m_N( AA_CACB.GetAtom( GetAtomTypes().N)),
      m_CA( AA_CACB.GetCA()),
      m_C( AA_CACB.GetAtom( GetAtomTypes().C)),
      m_O( AA_CACB.GetAtom( GetAtomTypes().O)),
      m_FirstSidechainAtom( AA_CACB.GetFirstSidechainAtom()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_N, m_CA, m_C, m_O, m_FirstSidechainAtom))
    {
    }

    //! @brief virtual copy constructor
    AABackBone *AABackBone::Clone() const
    {
      return new AABackBone( *this);
    }

    //! @brief virtual empty constructor with AAType
    //! @param SP_AA_DATA AAData object with information for the new AA
    //! this function is designed to be used in cases where AAClass is used for production multiple AA's
    //! which do not share one common AAData
    AABackBone *AABackBone::Empty( const util::ShPtr< AAData> &SP_AA_DATA) const
    {
      return new AABackBone( SP_AA_DATA);
    }

    //! @brief destructor
    AABackBone::~AABackBone()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AABackBone::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get types of atoms
    //! @return set of AtomTypes
    const storage::Set< AtomType> &AABackBone::GetTypesOfAtoms() const
    {
      // initialize a static set of AtomType from all backbone atoms
      static storage::Set< AtomType> s_atom_type_set
      (
        storage::Set< AtomType>::Create
        (
          GetAtomTypes().N,
          GetAtomTypes().CA,
          GetAtomTypes().C,
          GetAtomTypes().O,
          GetAtomTypes().CB,
          GetAtomTypes().HA2
        )
      );

      // return
      return s_atom_type_set;
    }

    //! @brief get all atoms of specified types
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return SiPtrVector of Atoms with specified type
    util::SiPtrVector< const Atom> AABackBone::GetAtoms( const storage::Set< AtomType> &ATOM_TYPES) const
    {
      // create and initialize a siptr vector
      util::SiPtrVector< const Atom> atoms;

      // allocate memory
      atoms.AllocateMemory( ATOM_TYPES.GetSize());

      // iterate over atom types
      for
      (
        storage::Set< AtomType>::const_iterator atom_itr( ATOM_TYPES.Begin()), atom_itr_end( ATOM_TYPES.End());
        atom_itr != atom_itr_end; ++atom_itr
      )
      {
        // get the atom
        const Atom &this_atom( GetAtom( *atom_itr));

        // only if the atom is defined
        if( this_atom.GetType().IsDefined())
        {
          // pushback the atom
          atoms.PushBack( util::ToSiPtr( this_atom));
        }
      }
      // return coordinates
      return atoms;
    }

    //! @brief set the specified atom
    //! @param ATOM Atom to be set
    void AABackBone::SetAtom( const Atom &ATOM)
    {
      // switch over type of ATOM and assign it to corresponding backbone atom with same type
           if( ATOM.GetType() == GetAtomTypes().N)  { m_N  = ATOM;}
      else if( ATOM.GetType() == GetAtomTypes().CA) { m_CA = ATOM;}
      else if( ATOM.GetType() == GetAtomTypes().C)  { m_C  = ATOM;}
      else if( ATOM.GetType() == GetAtomTypes().O)  { m_O  = ATOM;}
      else if( ATOM.GetType() == GetType()->GetFirstSidechainAtomType()) { m_FirstSidechainAtom = ATOM;}
      else
      {
        // if no match was found issue warning
        BCL_MessageStd
        (
          "This amino acid does not have the specified atom type " + ATOM.GetType().GetName()
        );
      }
    }

    //! @brief set all atoms
    //! @param ATOMS SiPtrVector of Atoms to be set
    void AABackBone::SetAtoms( const util::SiPtrVector< const Atom> &ATOMS)
    {
      // search for each individual backbone atom, and if it is found assign
      const util::SiPtr< const Atom> sp_n( Atom::FindAtom( ATOMS, GetAtomTypes().N));
      if( sp_n->GetType().IsDefined())
      {
        m_N = *sp_n;
      }
      const util::SiPtr< const Atom> sp_ca( Atom::FindAtom( ATOMS, GetAtomTypes().CA));
      if( sp_ca->GetType().IsDefined())
      {
        m_CA = *sp_ca;
      }
      const util::SiPtr< const Atom> sp_c( Atom::FindAtom( ATOMS, GetAtomTypes().C));
      if( sp_c->GetType().IsDefined())
      {
        m_C = *sp_c;
      }
      const util::SiPtr< const Atom> sp_o( Atom::FindAtom( ATOMS, GetAtomTypes().O));
      if( sp_o->GetType().IsDefined())
      {
        m_O = *sp_o;
      }
      const util::SiPtr< const Atom> sp_first_sc_atom( Atom::FindAtom( ATOMS, GetType()->GetFirstSidechainAtomType()));
      if( sp_first_sc_atom->GetType().IsDefined())
      {
        m_FirstSidechainAtom = *sp_first_sc_atom;
      }
    }

    //! @brief get AAClass
    //! @return AAClass
    const AAClass &AABackBone::GetAAClass() const
    {
      return GetAAClasses().e_AABackBone;
    }

    //! @brief calculate Omega backbone angle
    //! @param PREVIOUS_CA previous CA atom
    //! @param PREVIOUS_C previous carbon atom
    //! @return omega angle
    double AABackBone::CalculateOmega( const Atom &PREVIOUS_CA, const Atom &PREVIOUS_C) const
    {
      // if the type of PREVIOUS_CA is not CA
      if( PREVIOUS_CA.GetType() != GetAtomTypes().CA)
      {
        // issue warning
        BCL_MessageCrt
        (
          "given PREVIOUS_CA is not of type CA. but: " + util::Format()( PREVIOUS_CA.GetType())
        );

        // return undefined value
        return util::GetUndefinedDouble();
      }

      // if PREVIOUS_C is not a carbon
      if( PREVIOUS_C.GetType() != GetAtomTypes().C)
      {
        // issue warning
        BCL_MessageCrt
        (
          "given PREVIOUS_C is not of type C. but: " + util::Format()( PREVIOUS_C.GetType())
        );

        // return undefined
        return util::GetUndefinedDouble();
      }

      // calculate the dihedral angle and return it
      return Dihedral( PREVIOUS_CA, PREVIOUS_C, m_N, m_CA);
    }

    //! @brief calculate phi backbone angle if hydrogen is part of this aa
    //! @return phi angle
    double AABackBone::Phi() const
    {
      return util::GetUndefined< double>();
    }

    //! @brief calculate phi backbone angle
    //! @param PREVIOUS_C previous carbon atom
    //! @return phi angle
    double AABackBone::CalculatePhi( const Atom &PREVIOUS_C) const
    {
      // if PREVIOUS_C is not of type Carbon
      if( PREVIOUS_C.GetType() != GetAtomTypes().C)
      {
        // issue warning
        BCL_MessageCrt
        (
          "given PREVIOUS_C is not of type C. but: " + util::Format()( PREVIOUS_C.GetType())
        );

        // return undefined
        return util::GetUndefinedDouble();
      }

      // calculate dihedral and return it
      return Dihedral( PREVIOUS_C, m_N, m_CA, m_C);
    }

    //! @brief calculate psi backbone angle
    //! @param FOLLOWING_N following nitrogen atom
    //! @return psi angle
    double AABackBone::CalculatePsi( const Atom &FOLLOWING_N) const
    {
      // if FOLLOWING_N is not of type Nitrogen
      if( FOLLOWING_N.GetType() != GetAtomTypes().N)
      {
        // issue warning
        BCL_MessageCrt
        (
          "given FOLLOWING_N is not of type N. but: " + util::Format()( FOLLOWING_N.GetType())
        );

        // return undefined
        return util::GetUndefinedDouble();
      }

      // calculate dihedral and return it
      return Dihedral( m_N, m_CA, m_C, FOLLOWING_N);
    }

    //! @brief calculate psi backbone angle if oxygen is part of this aa
    //! @return psi angle
    double AABackBone::Psi() const
    {
      double psi( Dihedral( m_N, m_CA, m_C, m_O));

      // rotate by 180, since the oxygen is opposite to the nitrogen usually used for psi
      if( psi > 0.0)
      {
        psi -= math::g_Pi;
      }
      else
      {
        psi += math::g_Pi;
      }

      return psi;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all atom coordinates
    //! @return all atom coordinates
    util::SiPtrVector< const linal::Vector3D> AABackBone::GetAtomCoordinates() const
    {
      // create and initialize coordinates vector
      util::SiPtrVector< const linal::Vector3D> coordinates;
      coordinates.PushBack( util::ToSiPtr( m_N.GetCoordinates()));
      coordinates.PushBack( util::ToSiPtr( m_CA.GetCoordinates()));
      coordinates.PushBack( util::ToSiPtr( m_C.GetCoordinates()));
      coordinates.PushBack( util::ToSiPtr( m_O.GetCoordinates()));
      coordinates.PushBack( util::ToSiPtr( m_FirstSidechainAtom.GetCoordinates()));

      // return coordinates
      return coordinates;
    }

    //! @brief get all atom coordinates for the specified atom types
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return coordinates of atoms of types specified in ATOM_TYPES
    util::SiPtrVector< const linal::Vector3D> AABackBone::GetAtomCoordinates
    (
      const storage::Set< AtomType> &ATOM_TYPES
    ) const
    {
      // create and initialize a siptr vector
      util::SiPtrVector< const linal::Vector3D> coordinates;

      // allocate memory
      coordinates.AllocateMemory( ATOM_TYPES.GetSize());

      // iterate over atom types
      for
      (
        storage::Set< AtomType>::const_iterator atom_itr( ATOM_TYPES.Begin()), atom_itr_end( ATOM_TYPES.End());
        atom_itr != atom_itr_end; ++atom_itr
      )
      {
        // get the atom
        const Atom &this_atom( GetAtom( *atom_itr));

        // only if the atom is defined
        if( this_atom.GetType().IsDefined())
        {
          // pushback the coordinates
          coordinates.PushBack( util::ToSiPtr( this_atom.GetCoordinates()));
        }
      }
      // return coordinates
      return coordinates;
    }

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION translation vector to be applied
    void AABackBone::Translate( const linal::Vector3D &TRANSLATION)
    {
      // transform all atoms
      m_N.Translate( TRANSLATION);
      m_CA.Translate( TRANSLATION);
      m_C.Translate( TRANSLATION);
      m_O.Translate( TRANSLATION);
      m_FirstSidechainAtom.Translate( TRANSLATION);
    }

    //! @brief transform the object by a given TRANSFORMATION_MATRIX_3D
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void AABackBone::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      // transform all atoms
      m_N.Transform( TRANSFORMATION_MATRIX_3D);
      m_CA.Transform( TRANSFORMATION_MATRIX_3D);
      m_C.Transform( TRANSFORMATION_MATRIX_3D);
      m_O.Transform( TRANSFORMATION_MATRIX_3D);
      m_FirstSidechainAtom.Transform( TRANSFORMATION_MATRIX_3D);
    }

    //! @brief rotate the object by a given ROTATION_MATRIX_3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
    void AABackBone::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      // transform all atoms
      m_N.Rotate( ROTATION_MATRIX_3D);
      m_CA.Rotate( ROTATION_MATRIX_3D);
      m_C.Rotate( ROTATION_MATRIX_3D);
      m_O.Rotate( ROTATION_MATRIX_3D);
      m_FirstSidechainAtom.Rotate( ROTATION_MATRIX_3D);
    }

    //! @brief returns the geometric center of the object
    //! @return geometric center of the object
    linal::Vector3D AABackBone::GetCenter() const
    {
      return coord::CenterOfMass( GetAtomCoordinates());
    }

    //! @brief set the aminoacid to the ideal conformation according to the SS_TYPE with given TRANSFORMATION_MATRIX_3D
    //! @param SS_TYPE SSType this AABase derived class is in
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void AABackBone::SetToIdealConformation
    (
      const SSType &SS_TYPE,
      const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D
    )
    {
      // if SS_TYPE is not helix nor strand, then skip
      if( !SS_TYPE->IsStructured())
      {
        return;
      }

      // set each atom to corresponding SS_TYPE specific coordinates
      m_N.SetCoordinates
      (
        linal::Vector3D
        (
          GetAtomTypes().N->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
        ).Transform( TRANSFORMATION_MATRIX_3D)
      );
      m_CA.SetCoordinates
      (
        linal::Vector3D
        (
          GetAtomTypes().CA->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
        ).Transform( TRANSFORMATION_MATRIX_3D)
      );
      m_C.SetCoordinates
      (
        linal::Vector3D
        (
          GetAtomTypes().C->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
        ).Transform( TRANSFORMATION_MATRIX_3D)
      );
      m_O.SetCoordinates
      (
        linal::Vector3D
        (
          GetAtomTypes().O->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
        ).Transform( TRANSFORMATION_MATRIX_3D)
      );
      m_FirstSidechainAtom.SetCoordinates
      (
        linal::Vector3D
        (
          GetType()->GetFirstSidechainAtomType()->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
        ).Transform( TRANSFORMATION_MATRIX_3D)
      );
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator for AABackBone class
    //! @param AA_BACKBONE_RHS AABackBone object to be copied
    //! @return this after assignment to AA_BACKBONE_RHS is done
    AABackBone &AABackBone::operator =( const AABackBone &AA_BACKBONE_RHS)
    {
      // assign base class AABase and data members
      AABase::operator =( AA_BACKBONE_RHS);
      m_N = AA_BACKBONE_RHS.m_N;
      m_CA = AA_BACKBONE_RHS.m_CA;
      m_C = AA_BACKBONE_RHS.m_C;
      m_O = AA_BACKBONE_RHS.m_O;
      m_FirstSidechainAtom = AA_BACKBONE_RHS.m_FirstSidechainAtom;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AABackBone::Read( std::istream &ISTREAM)
    {
      // read base class
      AABase::Read( ISTREAM);

      // read members
      io::Serialize::Read( m_N, ISTREAM);
      io::Serialize::Read( m_CA, ISTREAM);
      io::Serialize::Read( m_C, ISTREAM);
      io::Serialize::Read( m_O, ISTREAM);
      io::Serialize::Read( m_FirstSidechainAtom, ISTREAM);

      //return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &AABackBone::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write base class
      AABase::Write( OSTREAM, INDENT) << '\n';

      // write members
      io::Serialize::Write( m_N, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_C, OSTREAM, INDENT)  << '\n';
      io::Serialize::Write( m_O, OSTREAM, INDENT)  << '\n';
      io::Serialize::Write( m_FirstSidechainAtom, OSTREAM, INDENT);

      //return
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_aa_base.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "align/bcl_align.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "biol/bcl_biol_aa.h"
#include "biol/bcl_biol_aa_classes.h"
#include "biol/bcl_biol_rotamer.h"
#include "chemistry/bcl_chemistry_atom_types.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AABase::AABase() :
      m_Data( util::ShPtr< AAData>( new AAData()))
    {
    }

    //! @brief construct AABase from util::ShPtr to AAData
    //! @brief SP_AA_DATA ShPTr to AAData to be copied
    AABase::AABase( const util::ShPtr< AAData> &SP_AA_DATA) :
      m_Data( SP_AA_DATA)
    {
    }

    //! @brief copy constructor - makes just a soft copy of m_Data
    //! @param AA_BASE_RHS AABase to be copied
    AABase::AABase( const AABase &AA_BASE_RHS) :
      m_Data( AA_BASE_RHS.m_Data)
    {
    }

    //! destructor
    AABase::~AABase()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get all atoms of the specified element type
    //! @param ELEMENT_TYPE element types of interest
    //! @return all atoms of the specified element type
    util::SiPtrVector< const Atom> AABase::GetAtoms( const chemistry::ElementType &ELEMENT_TYPE) const
    {
      // get all the atoms
      const util::SiPtrVector< const Atom> all_atoms( GetAtoms());
      util::SiPtrVector< const Atom> found_atoms;

      // iterate over the atoms
      for
      (
        util::SiPtrVector< const Atom>::const_iterator atom_itr( all_atoms.Begin()), atom_itr_end( all_atoms.End());
        atom_itr != atom_itr_end; ++atom_itr
      )
      {
        // if the atom is the right element type
        if( ( *atom_itr)->GetType()->GetElementType() == ELEMENT_TYPE)
        {
          // add it to the vector
          found_atoms.PushBack( *atom_itr);
        }
      }

      // end
      return found_atoms;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the center of mass for all side chain atoms
    //! @return center of mass (not weighted by mass) of all side chain atoms - undefined if there are no atoms
    linal::Vector3D AABase::CalculateCenterOfMassOfSideChain() const
    {
       // if the aa passed wasn't an AAComplete, then just return the center of the first side chain atom
      if( GetAAClass() != GetAAClasses().e_AAComplete)
      {
        // center of the first side chain atom is center of mass for the sidechain
        return GetFirstSidechainAtom().GetCoordinates();
      }

      // all the side chain atom types of that particular amino acid
      std::vector< AtomType> atom_types_sidechain;

      // fill atom types with difference of all atom types for that amino acid and the backbone amino acid types
      std::set_difference
      (
        GetTypesOfAtoms().Begin(), GetTypesOfAtoms().End(),
        GetAtomTypes().GetBackBoneAtomTypes().Begin(), GetAtomTypes().GetBackBoneAtomTypes().End(),
        std::back_insert_iterator< std::vector< AtomType> >( atom_types_sidechain)
      );

      // obtain all the coordinates for the side chain atoms
      util::SiPtrVector< const linal::Vector3D> side_chain_atom_coordinates
      (
        GetAtomCoordinates( storage::Set< AtomType>( atom_types_sidechain.begin(), atom_types_sidechain.end()))
      );

      // calculate center of mass of the side chain based on all the side chain coordinates
      return coord::CenterOfMass( side_chain_atom_coordinates);
    }

    //! @brief calculate phi and psi
    //! @param PREVIOUS_C previous carbon atom
    //! @param FOLLOWING_N followin nitrogen atom
    //! @return phi and psi angles
    const storage::Pair< double, double> AABase::CalculatePhiPsi( const Atom &PREVIOUS_C, const Atom &FOLLOWING_N) const
    {
      return storage::Pair< double, double>( CalculatePhi( PREVIOUS_C), CalculatePsi( FOLLOWING_N));
    }

    //! @brief create a locator for that aa
    //! @param USE_PDB_ID
    //! @return locator aa to locate an amino acid in a protein model
    assemble::LocatorAA AABase::Locator( const bool USE_PDB_ID) const
    {
      return assemble::LocatorAA( m_Data->GetChainID(), USE_PDB_ID ? m_Data->GetPdbID() : m_Data->GetSeqID(), USE_PDB_ID);
    }

  ////////////////////////////////////////
  // data access - SSPrediction related //
  ////////////////////////////////////////

    //! @brief return map of secondary structure predictions stored for various method
    //! @return map of secondary structure predictions stored for various method
    const storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> > &AABase::GetSSPredictions() const
    {
      return m_Data->m_SSPredictions;
    }

    //! @brief secondary structure predictions stored for given SS_METHOD
    //! @brief METHOD sspred::Method of interest
    //! @return secondary structure predictions stored for SS_METHOD
    util::SiPtr< const sspred::MethodInterface> AABase::GetSSPrediction
    (
      const sspred::Method &SS_METHOD
    ) const
    {
      return m_Data->GetSSPrediction( SS_METHOD);
    }

    //! @brief provide sspredictions from different methods in PREDICTIONS_MAP to overwrite existing ones
    //! @param PREDICTIONS_MAP Map of predictions with different methods
    void AABase::SetSSPredictions
    (
      const storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> > &PREDICTIONS_MAP
    )
    {
      m_Data->m_SSPredictions = PREDICTIONS_MAP;
    }

    //! @brief Sets the predictions for the given SS_METHOD
    //! @param SS_METHOD SSMethod of interest
    //! @param PREDICTION MethodInterface derived class that stores the predictions
    void AABase::SetSSPrediction
    (
      const sspred::Method &SS_METHOD,
      const sspred::MethodInterface &PREDICTION
    )
    {
      m_Data->m_SSPredictions[ SS_METHOD] = util::ShPtr< sspred::MethodInterface>( PREDICTION.Clone());
    }

    //! @brief remove structure-based secondary-structure and TM methods from the map for this AA
    void AABase::RemoveStructureBasedSSTMInfo()
    {
      for
      (
        storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> >::iterator
          itr( m_Data->m_SSPredictions.Begin()), itr_end( m_Data->m_SSPredictions.End());
        itr != itr_end;
        ++itr
      )
      {
        if( ( *itr->first)->GetIsDeterminedFromSturcture())
        {
          itr->second.Reset();
        }
      }
    }

    //! @brief get the exposure prediction
    //! @return exposure prediction
    const double &AABase::GetExposurePrediction() const
    {
      return m_Data->m_ExposurePrediction;
    }

    //! @brief set the exposure prediction
    //! @param EXPOSURE exposure to set
    void AABase::SetExposurePrediction( const double &EXPOSURE)
    {
      m_Data->m_ExposurePrediction = EXPOSURE;
    }

  /////////////////////////////////
  // data access - Blast related //
  /////////////////////////////////

    //! @brief returns BlastProfile Object
    //! @param ISTREAM input stream to read the blast profile from
    void AABase::ReadBlastProfile( std::istream &ISTREAM)
    {
      if( !m_Data->m_BlastProfile.IsDefined())
      {
        m_Data->m_BlastProfile = util::ShPtr< BlastProfile>( new BlastProfile);
      }
      m_Data->m_BlastProfile->ReadProfile( ISTREAM);
    }

    //! @brief returns const BlastProfile Object
    //! @return const BlastProfile Object
    const BlastProfile &AABase::GetBlastProfile() const
    {
      static const BlastProfile s_undefined_profile;
      return m_Data->m_BlastProfile.IsDefined() ? *m_Data->m_BlastProfile : s_undefined_profile;
    }

    //! @brief returns const BlastProfile Object
    //! @return const BlastProfile Object
    const util::ShPtr< BlastProfile> &AABase::GetBlastProfilePtr() const
    {
      return m_Data->m_BlastProfile;
    }

    //! @brief set BlastProfile and Probabilities accordingly to given BLAST_PROFILE object
    //! @param BLAST_PROFILE blast profile object to be set to
    void AABase::SetBlastProfile( const BlastProfile &BLAST_PROFILE)
    {
      m_Data->m_BlastProfile = util::ShPtr< BlastProfile>( BLAST_PROFILE.Clone());
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns whether all coordinates for the given atom types are defined
    //! @param ATOM_TYPES Atom Types of interest
    bool AABase::HasDefinedCoordinates( const storage::Set< AtomType> &ATOM_TYPES) const
    {
      // get all the atom coordinates
      const util::SiPtrVector< const linal::Vector3D> coordinates( GetAtomCoordinates( ATOM_TYPES));

      // iterate over all coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr( coordinates.Begin()), coord_itr_end( coordinates.End());
        coord_itr != coord_itr_end; ++coord_itr
      )
      {
        // if any one of them is undefined
        if( !( *coord_itr)->IsDefined())
        {
          // then return failure
          return false;
        }
      }

      // this point is reached only if all coordinates were defined, therefore return true
      return true;
    }

    //! @brief check if this aa precedes a given amino acid
    //! @param AMINO_ACID second amino acid, that could precede this
    //! @return true if both amino acids have same chain id and if (other seq id) - (this seq_id) == 1
    bool AABase::DoesPrecede( const AABase &AMINO_ACID) const
    {
      return AMINO_ACID.GetChainID() == GetChainID() && ( AMINO_ACID.GetSeqID() - GetSeqID()) == 1;
    }

    //! @brief gives dihedral angles of the side chain - starting closest to backbone and moving out along side chain
    //! @return vector holding the dihedral angles - first element is dihedral angle closest to backbone
    Rotamer AABase::CalculateSideChainDihedralAngles() const
    {
      // get the atom types that are involved in the side chain dihedral angles for this residue
      const storage::Map< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> > atom_types
      (
        m_Data->GetType()->GetSideChainDihedralAngleAtomTypes()
      );

      // will hold the dihedral angles
      Rotamer dihedral_angles;

      // true if there are not enough atom types to calculate a dihedral angle
      if( atom_types.IsEmpty())
      {
        BCL_MessageCrt
        (
          "cannot calculate dihedral angles of residue " + GetIdentification() +
          " because the atom types specified for calculating dihedral angles are " + util::Format()( atom_types)
        );
        // return empty vector
        return dihedral_angles;
      }

      // iterate through the atom types
      for
      (
        storage::Map< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >::const_iterator
          chi_itr( atom_types.Begin()), chi_itr_end( atom_types.End());
        chi_itr != chi_itr_end;
        ++chi_itr
      )
      {
        // get the atoms
        const linal::Vector3D &coords_a( GetAtom( chi_itr->second( 0)).GetCoordinates());
        const linal::Vector3D &coords_b( GetAtom( chi_itr->second( 1)).GetCoordinates());
        const linal::Vector3D &coords_c( GetAtom( chi_itr->second( 2)).GetCoordinates());
        const linal::Vector3D &coords_d( GetAtom( chi_itr->second( 3)).GetCoordinates());

        // true if any of the coordinates are not defined
        if( !coords_a.IsDefined() || !coords_b.IsDefined() || !coords_c.IsDefined() || !coords_d.IsDefined())
        {
          // add the dihedral to vector
          dihedral_angles.Insert( ChiAngle( chi_itr->first, util::GetUndefinedDouble(), math::Angle::e_Radian));

          // go to to next dihedral
          continue;
        }

        // get the dihedral angle
        const double dihedral( linal::Dihedral( coords_a, coords_b, coords_c, coords_d));

        // add the dihedral to vector
        dihedral_angles.Insert( ChiAngle( chi_itr->first, dihedral, math::Angle::e_Radian));
      }

      return dihedral_angles;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator for AABase class
    //! @param AA_BASE_RHS AABase to be assigned from
    //! @return this AABase after being updated to AA_BASE_RHS
    AABase &AABase::operator =( const AABase &AA_BASE_RHS)
    {
      // update data
      m_Data = AA_BASE_RHS.m_Data;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read fasta from given ISTREAM with given SEQ_ID
    //! @param ISTREAM input stream
    //! @param SEQ_ID seq id for this AABase
    //! @return std::istream from which was read
    std::istream &AABase::ReadFasta( std::istream &ISTREAM, const size_t SEQ_ID)
    {
      // initialize one_letter_code
      char one_letter_code;

      // skip newline and or empty characters
      do
      {
        // if the end of stream is reached
        if( !( ISTREAM >> one_letter_code))
        {
          // end
          return ISTREAM;
        }
      }
      while( one_letter_code == '\n' || one_letter_code == ' ');

      // change this aa to the AAType that come from the one letter code and also update the seqids.
      m_Data = util::ShPtr< AAData>
      (
        new AAData
        (
          GetAATypes().AATypeFromOneLetterCode( one_letter_code),
          SEQ_ID,
          SEQ_ID,
          m_Data->m_PdbICode,
          m_Data->m_ChainID
        )
      );

      //end
      return ISTREAM;
    }

    //! @brief write fasta to provided OSTREAM
    //! @param OSTREAM output stream
    //! @return std::ostream to which was written
    std::ostream &AABase::WriteFasta( std::ostream &OSTREAM) const
    {
      // write one letter code to OSTREAM
      OSTREAM << m_Data->GetType()->GetOneLetterCode();

      // return
      return OSTREAM;
    }

    //! @brief get Identification of this amino acid
    //! @return string with identification
    std::string AABase::GetIdentification() const
    {
      // initialize identification with sequence id, one letter and three letter code and one state ss prediction

      std::string identification
      (
        util::Format().W( 5)( GetSeqID()) + " " +
        GetType()->GetOneLetterCode() + " " +
        GetType()->GetThreeLetterCode() + " "
      );

      // check if the PDB SSType was set
      if( GetSSPrediction( sspred::GetMethods().e_PDB).IsDefined())
      {
        identification += GetSSPrediction( sspred::GetMethods().e_PDB)->GetOneStateSSPrediction()->GetOneLetterCode();
      }
      // if jufo was provided
      else if( GetSSPrediction( sspred::GetMethods().e_JUFO).IsDefined())
      {
        identification += GetSSPrediction( sspred::GetMethods().e_JUFO)->GetOneStateSSPrediction()->GetOneLetterCode();
      }
      // else print
      else
      {
        identification += "U";
      }

      // end
      return identification;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AABase::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Data, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write Identification to given OSTREAM
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return std::ostream to which was written
    std::ostream &AABase::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      //end
      return OSTREAM;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief Calculate the distance between two aminoacids by calling AtomDistance function for the first side chain atoms
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @return CB distance between AMINO_ACID_A and AMINO_ACID_B
    double FirstSidechainAtomDistance( const AABase &AMINO_ACID_A, const AABase &AMINO_ACID_B)
    {
      return Distance( AMINO_ACID_A.GetFirstSidechainAtom(), AMINO_ACID_B.GetFirstSidechainAtom());
    }

    //! @brief Calculate the sequence separation between two amino acids
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @return sequence separation between AMINO_ACID_A and AMINO_ACID_B, will be undefined if they are from the same
    //!          chain or if the sequence ids are the same
    size_t SequenceSeparation( const AABase &AMINO_ACID_A, const AABase &AMINO_ACID_B)
    {
      // if the chain ids are different return undefined
      if( AMINO_ACID_A.GetChainID() != AMINO_ACID_B.GetChainID())
      {
        return util::GetUndefined< size_t>();
      }
      // if the sequence ids are the same return undefined
      if( AMINO_ACID_A.GetSeqID() == AMINO_ACID_B.GetSeqID())
      {
        return util::GetUndefined< size_t>();
      }
      // if AMINO_ACID_A comes first
      if( AMINO_ACID_A.GetSeqID() < AMINO_ACID_B.GetSeqID())
      {
        return AMINO_ACID_B.GetSeqID() - AMINO_ACID_A.GetSeqID() - 1;
      }
      // else
      return AMINO_ACID_A.GetSeqID() - AMINO_ACID_B.GetSeqID() - 1;
    }

    //! @brief check if two amino acids are in range to be peptide bonded
    //! @param AA_LEFT left in sequence
    //! @param AA_RIGHT right in sequence
    //! @param TEST_OMEGA check the angle also
    //! @return true if amino acids are close enough for a peptide bond
    bool AABase::AreAminoAcidsPeptideBonded
    (
      const AABase &AA_LEFT,
      const AABase &AA_RIGHT,
      const bool TEST_OMEGA
    )
    {
      static const double s_peptide_bond_length( GetAtomTypes().C->GetBondLength( GetAtomTypes().N));
      static const double s_peptide_bond_tolerance( 0.02 * s_peptide_bond_length);

      const storage::VectorND< 2, double> bond_length_angle( PeptideBondLengthAndAngle( AA_LEFT, AA_RIGHT));

      // check the distance
      if( !util::IsDefined( bond_length_angle.First()) || !math::EqualWithinAbsoluteTolerance( s_peptide_bond_length, bond_length_angle.First(), s_peptide_bond_tolerance))
      {
//        if( bond_length_angle.First() < 1.5)
//        {
//          BCL_Message
//          (
//            util::Message::e_Verbose,
//            "barely missed peptide bond: " + util::Format()( bond_length_angle.First()) + " "
//            + AA_LEFT.GetIdentification() + " " + AA_RIGHT.GetIdentification()
//          );
//        }
        return false;
      }

      // distance matched, need to test angle also?
      if( !TEST_OMEGA)
      {
        return true;
      }

      // tolerance of 10 degrees
      static const double s_angle_tolerance( math::g_Pi / 18.0);

      // check angle
      if
      (
           !util::IsDefined( bond_length_angle.Second())                                             // will be undefined if one of the coordinates is not defined
        ||
           (
                 !math::EqualWithinAbsoluteTolerance(    -math::g_Pi, bond_length_angle.Second(), s_angle_tolerance) // check for cis peptide bond
              && !math::EqualWithinAbsoluteTolerance(            0.0, bond_length_angle.Second(), s_angle_tolerance) // check for cis peptide bond
              && !math::EqualWithinAbsoluteTolerance(     math::g_Pi, bond_length_angle.Second(), s_angle_tolerance) // check for trans peptide bond
           )
      )
      {
//        BCL_Message
//        (
//          util::Message::e_Verbose,
//          "missed peptide bond angle: " + util::Format()( bond_length_angle.Second()) + " " +
//          AA_LEFT.GetIdentification() + " " + AA_RIGHT.GetIdentification()
//        );
        return false;
      }

      // meets all conditions
      return true;
    }

    //! @brief calculate peptide bond length and angle
    //! @param AA_LEFT left in sequence
    //! @param AA_RIGHT right in sequence
    //! @return vector of peptide bond length (C->N) and peptide bond angle ca->c->n->ca)
    storage::VectorND< 2, double> AABase::PeptideBondLengthAndAngle
    (
      const AABase &AA_LEFT,
      const AABase &AA_RIGHT
    )
    {
      const linal::Vector3D &coord_ca_left( AA_LEFT.GetAtom( GetAtomTypes().CA).GetCoordinates());
      const linal::Vector3D &coord_c( AA_LEFT.GetAtom( GetAtomTypes().C).GetCoordinates());
      const linal::Vector3D &coord_n( AA_RIGHT.GetAtom( GetAtomTypes().N).GetCoordinates());
      const linal::Vector3D &coord_ca_right( AA_RIGHT.GetAtom( GetAtomTypes().CA).GetCoordinates());

      return
        storage::VectorND< 2, double>
        (
          linal::Distance( coord_c, coord_n), // bond
          linal::Dihedral( coord_ca_left, coord_c, coord_n, coord_ca_right) // angle
        );
    }

    //! @brief Compute the nearest atom separation between two AAs (nearest of atom distances - VdW radii of atoms involved)
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @return the nearest atom separation between two AAs (nearest of atom distances - VdW radii of atoms involved)
    double NearestAtomVdWSphereSeparation( const AABase &AMINO_ACID_A, const AABase &AMINO_ACID_B, const bool &IGNORE_BB_BB)
    {
      double nearest_so_far( std::numeric_limits< double>::max());
      storage::Vector< double> vdw_radii_b;
      vdw_radii_b.AllocateMemory( AMINO_ACID_B.GetAtoms().GetSize());
      const AAType type_a( AMINO_ACID_A.GetType()), type_b( AMINO_ACID_B.GetType());
      for( auto itr_b( AMINO_ACID_B.GetAtoms().Begin()), itr_b_end( AMINO_ACID_B.GetAtoms().End()); itr_b != itr_b_end; ++itr_b)
      {
        vdw_radii_b.PushBack
        (
          type_b->GetChemistryAtomType( ( *itr_b)->GetType())->GetAtomTypeProperty( chemistry::AtomTypeData::e_VdWaalsRadiusCSD)
        );
      }
      for( auto itr_a( AMINO_ACID_A.GetAtoms().Begin()), itr_a_end( AMINO_ACID_A.GetAtoms().End()); itr_a != itr_a_end; ++itr_a)
      {
        if( !( *itr_a)->GetCoordinates().IsDefined())
        {
          continue;
        }
        if( IGNORE_BB_BB && ( *itr_a)->GetType()->IsBackBone() && ( *itr_a)->GetType() != GetAtomTypes().CA)
        {
          continue;
        }
        const double vdw_radius_a
        (
          type_a->GetChemistryAtomType( ( *itr_a)->GetType())->GetAtomTypeProperty( chemistry::AtomTypeData::e_VdWaalsRadiusCSD)
        );
        auto itr_vdw_radii_b( vdw_radii_b.Begin());
        for( auto itr_b( AMINO_ACID_B.GetAtoms().Begin()), itr_b_end( AMINO_ACID_B.GetAtoms().End()); itr_b != itr_b_end; ++itr_b, ++itr_vdw_radii_b)
        {
          if( !( *itr_b)->GetCoordinates().IsDefined())
          {
            continue;
          }
          if( IGNORE_BB_BB && ( *itr_b)->GetType()->IsBackBone() && ( *itr_b)->GetType() != GetAtomTypes().CA)
          {
            continue;
          }
          nearest_so_far
            = std::min
              (
                nearest_so_far,
                std::max
                (
                  0.0,
                  linal::Distance( ( *itr_b)->GetCoordinates(), ( *itr_a)->GetCoordinates())
                  - *itr_vdw_radii_b - vdw_radius_a
                )
              );
        }
      }
      return nearest_so_far;
    }

    //! @brief construct single letter code and seqid
    //! @param ONE_LETTER_CODE aa one letter code
    //! @param SEQ_ID sequence id
    //! @return pointer to new amino acid
    AABase *AABase::Construct( const char ONE_LETTER_CODE, const int SEQ_ID)
    {
      // create ShPtr to AAData and set aatype, seq_id, chain_id
      util::ShPtr< AAData> new_member_data
      (
        new AAData
        (
          GetAATypes().AATypeFromOneLetterCode( toupper( ONE_LETTER_CODE)),
          SEQ_ID,
          AAData::s_DefaultPdbID,
          AAData::s_DefaultPdbICode,
          'A'
        )
      );

      // create AA
      return new AA( new_member_data);
    }

  } // namespace biol

  namespace align
  {
  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief specialization of the template function t_Member objects of type AABase
    //! @param MEMBER the t_Member object to extract the character id from
    //! @return the character id
    template<>
    char GetCharId< biol::AABase>( const biol::AABase &MEMBER)
    {
      return MEMBER.GetType()->GetOneLetterCode();
    }

    //! @brief specialization of the template function t_Member objects of type AABase
    //! @param MEMBER the t_Member object to extract the complete identifier from
    //! @return the complete identifier
    template<>
    std::string GetCompleteId< biol::AABase>( const biol::AABase &MEMBER)
    {
      return MEMBER.GetIdentification();
    }

  } // namespace align
} // namespace bcl
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
#include "biol/bcl_biol_aa_ca_cb.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_classes.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AACaCb::s_Instance
    (
      GetObjectInstances().AddInstance( new AACaCb())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AACaCb::AACaCb() :
      AABase(),
      m_CA( GetAtomTypes().CA),
      m_FirstSidechainAtom(),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_CA, m_FirstSidechainAtom))
    {
    }

    //! @brief construct AACaCb from util::ShPtr to AAData
    //! @brief SP_AA_DATA ShPTr to AAData to be copied
    AACaCb::AACaCb( const util::ShPtr< AAData> &SP_AA_DATA) :
      AABase( SP_AA_DATA),
      m_CA( GetAtomTypes().CA),
      m_FirstSidechainAtom( SP_AA_DATA->GetType()->GetFirstSidechainAtomType()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_CA, m_FirstSidechainAtom))
    {
    }

    //! @brief constructor from AABase
    //! @param AA_BASE AABase
    AACaCb::AACaCb( const AABase &AA_BASE) :
      AABase( AA_BASE),
      m_CA( GetAtomTypes().CA),
      m_FirstSidechainAtom( GetType()->GetFirstSidechainAtomType()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_CA, m_FirstSidechainAtom))
    {
    }

    //! @brief construct AACaCb from AABase and atoms CA, CB
    //! @param AA_BASE AABase
    //! @param ATOM_CA CA atom
    //! @param ATOM_CB CB atom
    AACaCb::AACaCb
    (
      const AABase &AA_BASE,
      const Atom &ATOM_CA,
      const Atom &ATOM_CB
    ) :
      AABase( AA_BASE),
      m_CA( ATOM_CA),
      m_FirstSidechainAtom( ATOM_CB),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_CA, m_FirstSidechainAtom))
    {
    }

    //! @brief construct AACaCb from AABASE and util::SiPtrVector< const Atom> ATOMS
    //! @param AA_BASE AABase
    //! @param ATOMS SiPtrVector of Atoms
    AACaCb::AACaCb
    (
      const AABase &AA_BASE,
      const util::SiPtrVector< const Atom> &ATOMS
    ) :
      AABase( AA_BASE),
      m_CA( Atom::FindAtom( ATOMS, GetAtomTypes().CA)->GetType().IsDefined() ? *Atom::FindAtom( ATOMS, GetAtomTypes().CA) : GetAtomTypes().CA),
      m_FirstSidechainAtom( Atom::FindAtom( ATOMS, GetType()->GetFirstSidechainAtomType())->GetType().IsDefined() ? *Atom::FindAtom( ATOMS, GetType()->GetFirstSidechainAtomType()) : GetType()->GetFirstSidechainAtomType()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_CA, m_FirstSidechainAtom))
    {
    }

    //! @brief copy constructor
    //! @param AACACB AACaCb to be copied
    AACaCb::AACaCb( const AACaCb &AACACB) :
      AABase( AACACB),
      m_CA( AACACB.m_CA),
      m_FirstSidechainAtom( AACACB.m_FirstSidechainAtom),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_CA, m_FirstSidechainAtom))
    {
    }

    //! @brief virtual copy constructor
    AACaCb *AACaCb::Clone() const
    {
      return new AACaCb( *this);
    }

    //! @brief virtual empty constructor with AAData
    //! @param SP_AA_DATA AAData object with information for the new AA
    //! this function is designed to be used in cases where AAClass is used for production multiple AA's
    //! which do not share one common AAData
    AACaCb *AACaCb::Empty( const util::ShPtr< AAData> &SP_AA_DATA) const
    {
      return new AACaCb( SP_AA_DATA);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AACaCb::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get all atoms of specified types
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return SiPtrVector of Atoms with specified type
    util::SiPtrVector< const Atom> AACaCb::GetAtoms( const storage::Set< AtomType> &ATOM_TYPES) const
    {
      // create and initialize a siptr vector
      util::SiPtrVector< const Atom> atoms;

      // allocate memory
      atoms.AllocateMemory( ATOM_TYPES.GetSize());

      // iterate over atom types
      for
      (
        storage::Set< AtomType>::const_iterator atom_itr( ATOM_TYPES.Begin()), atom_itr_end( ATOM_TYPES.End());
        atom_itr != atom_itr_end; ++atom_itr
      )
      {
        // get the atom
        const Atom &this_atom( GetAtom( *atom_itr));

        // only if the atom is defined
        if( this_atom.GetType().IsDefined())
        {
          // pushback the atom
          atoms.PushBack( util::ToSiPtr( this_atom));
        }
      }
      // return coordinates
      return atoms;
    }

    //! @brief set the specified atom
    //! @param ATOM Atom to be set
    void AACaCb::SetAtom( const Atom &ATOM)
    {
      // if ATOM is CA set CA
      if( ATOM.GetType() == GetAtomTypes().CA)
      {
        m_CA = ATOM;
      }
      // else if ATOM is first side chain atom
      else if( ATOM.GetType() == GetType()->GetFirstSidechainAtomType())
      {
        m_FirstSidechainAtom = ATOM;
      }

      // if else issue warning
      BCL_MessageStd
      (
        "This amino acid does not have the specified atom type " + ATOM.GetType().GetName()
      );
    }

    //! @brief set all atoms
    //! @param ATOMS SiPtrVector of Atoms to be set
    void AACaCb::SetAtoms( const util::SiPtrVector< const Atom> &ATOMS)
    {
      // if CA is found assign CA
      const util::SiPtr< const Atom> sp_ca( Atom::FindAtom( ATOMS, GetAtomTypes().CA));
      if( sp_ca->GetType().IsDefined())
      {
        m_CA = *sp_ca;
      }
      // if first side chain atom is found, assign first side chain atom
      const util::SiPtr< const Atom> sp_cb( Atom::FindAtom( ATOMS, GetType()->GetFirstSidechainAtomType()));
      if( sp_cb->GetType().IsDefined())
      {
        m_FirstSidechainAtom = *sp_cb;
      }
    }

    //! @brief get AAClass
    //! @return AAClass
    const AAClass &AACaCb::GetAAClass() const
    {
      return GetAAClasses().e_AACaCb;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all atom coordinates
    //! @return all atom coordinates
    util::SiPtrVector< const linal::Vector3D> AACaCb::GetAtomCoordinates() const
    {

      // create and initialize coordinates vector
      util::SiPtrVector< const linal::Vector3D> coordinates;
      coordinates.PushBack( util::ToSiPtr( m_CA.GetCoordinates()));
      coordinates.PushBack( util::ToSiPtr( m_FirstSidechainAtom.GetCoordinates()));

      // return coordinates
      return coordinates;
    }

    //! @brief get all atom coordinates for the specified atom types
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return coordinates of atoms of types specified in ATOM_TYPES
    util::SiPtrVector< const linal::Vector3D> AACaCb::GetAtomCoordinates
    (
      const storage::Set< AtomType> &ATOM_TYPES
    ) const
    {
      // create and initialize a siptr vector
      util::SiPtrVector< const linal::Vector3D> coordinates;
      coordinates.AllocateMemory( ATOM_TYPES.GetSize());

      for
      (
        storage::Set< AtomType>::const_iterator atom_itr( ATOM_TYPES.Begin()), atom_itr_end( ATOM_TYPES.End());
        atom_itr != atom_itr_end; ++atom_itr
      )
      {
        // get the atom
        const Atom &this_atom( GetAtom( *atom_itr));

        // only if the atom is defined
        if( this_atom.GetType().IsDefined())
        {
          // pushback the coordinates
          coordinates.PushBack( util::ToSiPtr( this_atom.GetCoordinates()));
        }
      }
      // return coordinates
      return coordinates;
    }

    //! @brief set the aminoacid to the ideal conformation according to the SS_TYPE with given TRANSFORMATION_MATRIX_3D
    //! @param SS_TYPE SSType this AABase derived class is in
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void AACaCb::SetToIdealConformation
    (
      const SSType &SS_TYPE,
      const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D
    )
    {
      // if ss type is not helix or strand skip
      if( !SS_TYPE->IsStructured())
      {
        return;
      }

      // set CA coordinates to corresponding to SSType
      m_CA.SetCoordinates
      (
        linal::Vector3D
        (
          GetAtomTypes().CA->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
        ).Transform( TRANSFORMATION_MATRIX_3D)
      );

      // set first side chain atom coordinates to corresponding to SSType
      m_FirstSidechainAtom.SetCoordinates
      (
        linal::Vector3D
        (
          GetType()->GetFirstSidechainAtomType()->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
        ).Transform( TRANSFORMATION_MATRIX_3D)
      );
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator for AACaCb class
    //! @param AACACB_RHS AACaCb object to be copied
    //! @return this after assignment to AACACB_RHS is done
    AACaCb &AACaCb::operator =( const AACaCb &AACACB_RHS)
    {
      // assign base class and data members
      AABase::operator =( AACACB_RHS);
      m_CA = AACACB_RHS.m_CA;
      m_FirstSidechainAtom = AACACB_RHS.m_FirstSidechainAtom;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AACaCb::Read( std::istream &ISTREAM)
    {
      // read base class
      AABase::Read( ISTREAM);

      // read member
      io::Serialize::Read( m_CA, ISTREAM);
      io::Serialize::Read( m_FirstSidechainAtom, ISTREAM);

      //return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &AACaCb::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write base class
      AABase::Write( OSTREAM, INDENT) << '\n';

      // write members
      io::Serialize::Write( m_CA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FirstSidechainAtom, OSTREAM, INDENT);

      //return
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_aa_classes.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa.h"
#include "biol/bcl_biol_aa_back_bone.h"
#include "biol/bcl_biol_aa_complete.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct all AAClasses
    AAClasses::AAClasses() :
      e_AA(         AddEnum( "AA",         util::ShPtr< AABase>( new AA()))),
      e_AACaCb(     AddEnum( "AACaCb",     util::ShPtr< AABase>( new AACaCb()))),
      e_AABackBone( AddEnum( "AABackBone", util::ShPtr< AABase>( new AABackBone()))),
      e_AAComplete( AddEnum( "AAComplete", util::ShPtr< AABase>( new AAComplete())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAClasses::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief construct on access function for all AAClasses
    //! @return reference to only instances of AAClasses
    const AAClasses &GetAAClasses()
    {
      return AAClasses::GetEnums();
    }

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< biol::AABase>, biol::AAClasses>;

  } // namespace util
} // namespace bcl
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
#include "biol/bcl_biol_aa_compare.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////////////////
  // AACompareBySeqID //
  //////////////////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a seq id
    //! @param SEQ_ID Sequence id
    AACompareBySeqID::AACompareBySeqID( const int SEQ_ID) :
      m_SeqID( SEQ_ID)
    {
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief compares and returns whether the seq id of the given amino acid matched the stored one
    //! @param AMINO_ACID amino acoid to be compared
    //! @return returns whether the seq id of the given amino acid matched the stored one
    bool AACompareBySeqID::operator()( const AABase &AMINO_ACID) const
    {
      return AMINO_ACID.GetSeqID() == m_SeqID;
    }

    //! @brief compares and returns whether the seq id of the given amino acid matched the stored one
    //! @param SP_AMINO_ACID amino acoid to be compared
    //! @return returns whether the seq id of the given amino acid matched the stored one
    bool AACompareBySeqID::operator()( const util::PtrInterface< AABase> &SP_AMINO_ACID) const
    {
      return operator()( *SP_AMINO_ACID);
    }

    //! @brief compares and returns whether the seq id of the given amino acid matched the stored one
    //! @param SP_AMINO_ACID amino acoid to be compared
    //! @return returns whether the seq id of the given amino acid matched the stored one
    bool AACompareBySeqID::operator()( const util::PtrInterface< const AABase> &SP_AMINO_ACID) const
    {
      return operator()( *SP_AMINO_ACID);
    }

  //////////////////////
  // AACompareDataPtr //
  //////////////////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new class_name
    AACompareDataPtr *AACompareDataPtr::Clone() const
    {
      return new AACompareDataPtr( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AACompareDataPtr::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief () operator for returning whether two AABase are identical
    //! @param AA_BASE_LHS first AABase
    //! @param AA_BASE_RHS first AABase
    //! @return the whether two AABase are identical
    bool AACompareDataPtr::operator()( const AABase &AA_BASE_LHS, const AABase &AA_BASE_RHS) const
    {
      return AA_BASE_LHS.GetData() == AA_BASE_RHS.GetData();
    }

    //! @brief () operator for returning whether two AABase are identical
    //! @param SP_AA_BASE_LHS first AABase
    //! @param SP_AA_BASE_RHS first AABase
    //! @return the whether two AABase are identical
    bool AACompareDataPtr::operator()
    (
      const util::PtrInterface< const AABase> &SP_AA_BASE_LHS,
      const util::PtrInterface< const AABase> &SP_AA_BASE_RHS
    ) const
    {
      return operator ()( *SP_AA_BASE_LHS, *SP_AA_BASE_RHS);
    }

    //! @brief () operator for returning whether two AABase are identical
    //! @param SP_AA_BASE_LHS first AABase
    //! @param SP_AA_BASE_RHS first AABase
    //! @return the whether two AABase are identical
    bool AACompareDataPtr::operator()
    (
      const util::PtrInterface< AABase> &SP_AA_BASE_LHS,
      const util::PtrInterface< AABase> &SP_AA_BASE_RHS
    ) const
    {
      return operator ()( *SP_AA_BASE_LHS, *SP_AA_BASE_RHS);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AACompareDataPtr::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AACompareDataPtr::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  ///////////////////
  // AACompareData //
  ///////////////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new class_name
    AACompareData *AACompareData::Clone() const
    {
      return new AACompareData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AACompareData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief () operator for returning whether two AABase are identical
    //! @param AA_BASE_LHS first AABase
    //! @param AA_BASE_RHS first AABase
    //! @return the whether two AABase are identical
    bool AACompareData::operator()( const AABase &AA_BASE_LHS, const AABase &AA_BASE_RHS) const
    {
      return AA_BASE_LHS.GetChainID() == AA_BASE_RHS.GetChainID() &&
        AA_BASE_LHS.GetSeqID() == AA_BASE_RHS.GetSeqID() &&
        AA_BASE_LHS.GetType() == AA_BASE_RHS.GetType();
    }

    //! @brief () operator for returning whether two AABase are identical
    //! @param SP_AA_BASE_LHS first AABase
    //! @param SP_AA_BASE_RHS first AABase
    //! @return the whether two AABase are identical
    bool AACompareData::operator()
    (
      const util::PtrInterface< const AABase> &SP_AA_BASE_LHS,
      const util::PtrInterface< const AABase> &SP_AA_BASE_RHS
    ) const
    {
      return operator ()( *SP_AA_BASE_LHS, *SP_AA_BASE_RHS);
    }

    //! @brief () operator for returning whether two AABase are identical
    //! @param SP_AA_BASE_LHS first AABase
    //! @param SP_AA_BASE_RHS first AABase
    //! @return the whether two AABase are identical
    bool AACompareData::operator()
    (
      const util::PtrInterface< AABase> &SP_AA_BASE_LHS,
      const util::PtrInterface< AABase> &SP_AA_BASE_RHS
    ) const
    {
      return operator ()( *SP_AA_BASE_LHS, *SP_AA_BASE_RHS);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AACompareData::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AACompareData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  /////////////////////
  // AALessThanSeqID //
  /////////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief () operator for returning whether lhs amino acid is less than rhs amino acids by chain and seq id
    //! @param AA_BASE_LHS first AABase
    //! @param AA_BASE_RHS first AABase
    //! @return the whether one aa is less than the other aa
    bool AALessThanSeqID::operator()
    (
      const AABase &AA_BASE_LHS,
      const AABase &AA_BASE_RHS
    ) const
    {
      // compare chain id
      if( AA_BASE_LHS.GetChainID() == AA_BASE_RHS.GetChainID())
      {
        return AA_BASE_LHS.GetSeqID() < AA_BASE_RHS.GetSeqID();
      }

      return AA_BASE_LHS.GetChainID() < AA_BASE_RHS.GetChainID();
    }

    //! @brief () operator for returning whether lhs amino acid is less than rhs amino acids by chain and seq id
    //! @param SP_AA_BASE_LHS first AABase
    //! @param SP_AA_BASE_RHS first AABase
    //! @return the whether one aa is less than the other aa
    bool AALessThanSeqID::operator()
    (
      const util::PtrInterface< const AABase> &SP_AA_BASE_LHS,
      const util::PtrInterface< const AABase> &SP_AA_BASE_RHS
    ) const
    {
      // compare chain id
      if( SP_AA_BASE_LHS->GetChainID() == SP_AA_BASE_RHS->GetChainID())
      {
        return SP_AA_BASE_LHS->GetSeqID() < SP_AA_BASE_RHS->GetSeqID();
      }

      return SP_AA_BASE_LHS->GetChainID() < SP_AA_BASE_RHS->GetChainID();
    }

    //! @brief () operator for returning whether lhs amino acid is less than rhs amino acids by chain and seq id
    //! @param SP_AA_BASE_LHS first AABase
    //! @param SP_AA_BASE_RHS first AABase
    //! @return the whether one aa is less than the other aa
    bool AALessThanSeqID::operator()
    (
      const util::PtrInterface< AABase> &SP_AA_BASE_LHS,
      const util::PtrInterface< AABase> &SP_AA_BASE_RHS
    ) const
    {
      // compare chain id
      if( SP_AA_BASE_LHS->GetChainID() == SP_AA_BASE_RHS->GetChainID())
      {
        return SP_AA_BASE_LHS->GetSeqID() < SP_AA_BASE_RHS->GetSeqID();
      }

      return SP_AA_BASE_LHS->GetChainID() < SP_AA_BASE_RHS->GetChainID();
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_aa_complete.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_back_bone.h"
#include "biol/bcl_biol_aa_classes.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAComplete::s_Instance
    (
      GetObjectInstances().AddInstance( new AAComplete())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AAComplete::AAComplete() :
      AABase(),
      m_Atoms( DefaultBackboneAtoms()),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
    }

    //! @brief construct AABase from util::ShPtr to AAData
    //! @brief SP_AA_DATA ShPTr to AAData to be copied
    AAComplete::AAComplete( const util::ShPtr< AAData> &SP_AA_DATA) :
      AABase( SP_AA_DATA),
      m_Atoms( DefaultBackboneAtomsWithFirstSideChainAtom( SP_AA_DATA->GetType())),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
    }

    //! @brief constructor from AABase
    //! @param AA_BASE AABase
    AAComplete::AAComplete( const AABase &AA_BASE) :
      AABase( AA_BASE),
      m_Atoms( DefaultBackboneAtomsWithFirstSideChainAtom( AA_BASE.GetType())),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
      SetAtoms( AA_BASE.GetAtoms());
    }

    //! @brief construct AACaCb from AABASE and util::SiPtrVector< const Atom> ATOMS
    //! @param AA_BASE AABase
    //! @param ATOMS SiPtrVector of Atoms
    AAComplete::AAComplete( const AABase &AA_BASE, const util::SiPtrVector< const Atom> &ATOMS) :
      AABase( AA_BASE),
      m_Atoms( DefaultBackboneAtomsWithFirstSideChainAtom( AA_BASE.GetType())),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
      SetAtoms( ATOMS);
    }

    //! @brief copy constructor
    //! @param AA_COMPLETE AAComplete to be copied
    AAComplete::AAComplete( const AAComplete &AA_COMPLETE) :
      AABase( AA_COMPLETE),
      m_Atoms( DefaultBackboneAtomsWithFirstSideChainAtom( AA_COMPLETE.GetType())),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
      SetAtoms( AA_COMPLETE.GetAtoms());
    }

    //! @brief copy constructor
    //! @param AA_BACKBONE AABackBone to be copied
    AAComplete::AAComplete( const AABackBone &AA_BACKBONE) :
      AABase( AA_BACKBONE),
      m_Atoms( DefaultBackboneAtomsWithFirstSideChainAtom( AA_BACKBONE.GetType())),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
      SetAtoms( AA_BACKBONE.GetAtoms());
    }

    //! @brief constructor from CaCb amino acids
    //! @param AA_CACB AACaCb to be copied
    AAComplete::AAComplete( const AACaCb &AA_CACB) :
      AABase( AA_CACB),
      m_Atoms( DefaultBackboneAtomsWithFirstSideChainAtom( AA_CACB.GetType())),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
      SetAtoms( AA_CACB.GetAtoms());
    }

    //! @brief virtual copy constructor
    AAComplete *AAComplete::Clone() const
    {
      return new AAComplete( *this);
    }

    //! @brief virtual empty constructor with AAData
    //! @param SP_AA_DATA AAData object with information for the new AA
    //! this function is designed to be used in cases where AAClass is used for production multiple AA's
    //! which do not share one common AAData
    AAComplete *AAComplete::Empty( const util::ShPtr< AAData> &SP_AA_DATA) const
    {
      return new AAComplete( SP_AA_DATA);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAComplete::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get types of atoms
    //! @return set of AtomTypes
    const storage::Set< AtomType> &AAComplete::GetTypesOfAtoms() const
    {
      // return the atom types that belong to that complete amino acid
      return AABase::GetType()->GetAllowedAtomTypes();
    }

    //! @brief get all atoms of specified types
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return SiPtrVector of Atoms with specified type
    util::SiPtrVector< const Atom> AAComplete::GetAtoms( const storage::Set< AtomType> &ATOM_TYPES) const
    {
      // create and initialize a siptr vector
      util::SiPtrVector< const Atom> atoms;

      // allocate memory
      atoms.AllocateMemory( ATOM_TYPES.GetSize());

      // iterate over atom types
      for
      (
        storage::Set< AtomType>::const_iterator atom_itr( ATOM_TYPES.Begin()), atom_itr_end( ATOM_TYPES.End());
        atom_itr != atom_itr_end; ++atom_itr
      )
      {
        // get the atom
        const Atom &this_atom( GetAtom( *atom_itr));

        // only if the atom is defined
        if( this_atom.GetType().IsDefined())
        {
          // pushback the atom
          atoms.PushBack( util::ToSiPtr( this_atom));
        }
      }
      // return coordinates
      return atoms;
    }

    //! @brief get CB atom
    //! @return CB atom
    const Atom &AAComplete::GetFirstSidechainAtom() const
    {
      // return cb, it amino acid type has cb, HA2 otherwise (for GLY type amino acids)
      return GetAtom( GetType()->GetFirstSidechainAtomType());
    }

    //! @brief set the specified atom
    //! @param ATOM Atom to be set
    void AAComplete::SetAtom( const Atom &ATOM)
    {
      const AtomType &given_atom_type( ATOM.GetType());

      // iterate over all atoms to find the one with the correct type
      for
      (
        storage::List< Atom>::iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        // check if atom agrees
        if( itr->GetType() == given_atom_type)
        {
          *itr = ATOM;
          return;
        }
      }

      // if there was non found, check if the type does match the aa types atom types or is a terminal atom type and insert
      if( GetType()->DoesContainAtomType( given_atom_type) || GetAtomTypes().GetTerminalExtraAtomTypes().Has( given_atom_type))
      {
        m_Atoms.PushBack( ATOM);
      }
      // the atom cannot go into the amino acid
      else
      {
        // if no match was found issue warning
        BCL_MessageVrb
        (
          "This amino acid does not have the atom of the specified atom type " + ATOM.GetType().GetName() +
          "\nthe type of the amino acid type " + GetType().GetName() + " does not allow it"
        );

        // end
        return;
      }

      // update the siptr vector
      m_AtomList = util::SiPtrVector< const Atom>( m_Atoms.Begin(), m_Atoms.End());
    }

    //! @brief set all atoms
    //! @param ATOMS SiPtrVector of Atoms to be set
    void AAComplete::SetAtoms( const util::SiPtrVector< const Atom> &ATOMS)
    {
      // call set atom for each given atom
      for
      (
        util::SiPtrVector< const Atom>::const_iterator itr( ATOMS.Begin()), itr_end( ATOMS.End());
        itr != itr_end;
        ++itr
      )
      {
        SetAtom( **itr);
      }
    }

    //! @brief get AAClass
    //! @return AAClass
    const AAClass &AAComplete::GetAAClass() const
    {
      return GetAAClasses().e_AAComplete;
    }

    //! @brief calculate Omega backbone angle
    //! @param PREVIOUS_CA previous CA atom
    //! @param PREVIOUS_C previous carbon atom
    //! @return omega angle
    double AAComplete::CalculateOmega( const Atom &PREVIOUS_CA, const Atom &PREVIOUS_C) const
    {
      // if the type of PREVIOUS_CA is not CA
      if( PREVIOUS_CA.GetType() != GetAtomTypes().CA)
      {
        // issue warning
        BCL_MessageCrt
        (
          "given PREVIOUS_CA is not of type CA. but: " + util::Format()( PREVIOUS_CA.GetType())
        );

        // return undefined value
        return util::GetUndefinedDouble();
      }

      // if PREVIOUS_C is not a carbon
      if( PREVIOUS_C.GetType() != GetAtomTypes().C)
      {
        // issue warning
        BCL_MessageCrt
        (
          "given PREVIOUS_C is not of type C. but: " + util::Format()( PREVIOUS_C.GetType())
        );

        // return undefined
        return util::GetUndefinedDouble();
      }

      // calculate the dihedral angle and return it
      return Dihedral( PREVIOUS_CA, PREVIOUS_C, GetAtom( GetAtomTypes().N), GetAtom( GetAtomTypes().CA));
    }

    //! @brief calculate phi backbone angle if hydrogen is part of this aa
    //! @return phi angle
    double AAComplete::Phi() const
    {
      // calculate dihedral
      double phi( Dihedral( GetAtom( GetAtomTypes().H), GetAtom( GetAtomTypes().N), GetAtom( GetAtomTypes().CA), GetAtom( GetAtomTypes().C)));
      if( phi > 0.0)
      {
        phi -= math::g_Pi;
      }
      else
      {
        phi += math::g_Pi;
      }

      // return
      return phi;
    }

    //! @brief calculate phi backbone angle
    //! @param PREVIOUS_C previous carbon atom
    //! @return phi angle
    double AAComplete::CalculatePhi( const Atom &PREVIOUS_C) const
    {
      // if PREVIOUS_C is not of type Carbon
      if( PREVIOUS_C.GetType() != GetAtomTypes().C)
      {
        // issue warning
        BCL_MessageCrt
        (
          "given PREVIOUS_C is not of type C. but: " + util::Format()( PREVIOUS_C.GetType())
        );

        // return undefined
        return util::GetUndefinedDouble();
      }

      // calculate dihedral and return it
      return Dihedral( PREVIOUS_C, GetAtom( GetAtomTypes().N), GetAtom( GetAtomTypes().CA), GetAtom( GetAtomTypes().C));
    }

    //! @brief calculate psi backbone angle if oxygen is part of this aa
    //! @return psi angle
    double AAComplete::Psi() const
    {
      // calculate dihedral
      double psi( Dihedral( GetAtom( GetAtomTypes().N), GetAtom( GetAtomTypes().CA), GetAtom( GetAtomTypes().C), GetAtom( GetAtomTypes().O)));
      if( psi > 0.0)
      {
        psi -= math::g_Pi;
      }
      else
      {
        psi += math::g_Pi;
      }

      // return
      return psi;
    }

    //! @brief calculate psi backbone angle
    //! @param FOLLOWING_N following nitrogen atom
    //! @return psi angle
    double AAComplete::CalculatePsi( const Atom &FOLLOWING_N) const
    {
      // if FOLLOWING_N is not of type Nitrogen
      if( FOLLOWING_N.GetType() != GetAtomTypes().N)
      {
        // issue warning
        BCL_MessageCrt
        (
          "given FOLLOWING_N is not of type N. but: " + util::Format()( FOLLOWING_N.GetType())
        );

        // return undefined
        return util::GetUndefinedDouble();
      }

      // calculate dihedral and return it
      return Dihedral( GetAtom( GetAtomTypes().N), GetAtom( GetAtomTypes().CA), GetAtom( GetAtomTypes().C), FOLLOWING_N);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all atom coordinates
    //! @return all atom coordinates
    util::SiPtrVector< const linal::Vector3D> AAComplete::GetAtomCoordinates() const
    {
      // create and initialize coordinates vector
      util::SiPtrVector< const linal::Vector3D> coordinates;

      // iterate over all atoms and collect their coordinates
      for
      (
        storage::List< Atom>::const_iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        coordinates.PushBack( util::ToSiPtr( itr->GetCoordinates()));
      }

      // return coordinates
      return coordinates;
    }

    //! @brief get all atom coordinates for the specified atom types
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return coordinates of atoms of types specified in ATOM_TYPES
    util::SiPtrVector< const linal::Vector3D> AAComplete::GetAtomCoordinates
    (
      const storage::Set< AtomType> &ATOM_TYPES
    ) const
    {
      // create and initialize a siptr vector
      util::SiPtrVector< const linal::Vector3D> coordinates;

      // allocate memory
      coordinates.AllocateMemory( ATOM_TYPES.GetSize());

      // iterate over atom types
      for
      (
        storage::Set< AtomType>::const_iterator atom_itr( ATOM_TYPES.Begin()), atom_itr_end( ATOM_TYPES.End());
        atom_itr != atom_itr_end;
        ++atom_itr
      )
      {
        // get the atom
        const Atom &this_atom( GetAtom( *atom_itr));

        // only if the atom is defined
        if( this_atom.GetType().IsDefined())
        {
          // pushback the coordinates
          coordinates.PushBack( util::ToSiPtr( this_atom.GetCoordinates()));
        }
      }
      // return coordinates
      return coordinates;
    }

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION translation vector to be applied
    void AAComplete::Translate( const linal::Vector3D &TRANSLATION)
    {
      // transform all atoms
      for
      (
        storage::List< Atom>::iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->Translate( TRANSLATION);
      }
    }

    //! @brief transform the object by a given TRANSFORMATION_MATRIX_3D
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void AAComplete::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      // transform all atoms
      for
      (
        storage::List< Atom>::iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->Transform( TRANSFORMATION_MATRIX_3D);
      }
    }

    //! @brief rotate the object by a given ROTATION_MATRIX_3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
    void AAComplete::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      // transform all atoms
      for
      (
        storage::List< Atom>::iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->Rotate( ROTATION_MATRIX_3D);
      }
    }

    //! @brief returns the geometric center of the object
    //! @return geometric center of the object
    linal::Vector3D AAComplete::GetCenter() const
    {
      return coord::CenterOfMass( GetAtomCoordinates());
    }

    //! @brief set the aminoacid to the ideal conformation according to the SS_TYPE with given TRANSFORMATION_MATRIX_3D
    //! @param SS_TYPE SSType this AABase derived class is in
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void AAComplete::SetToIdealConformation
    (
      const SSType &SS_TYPE,
      const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D
    )
    {
      // if SS_TYPE is not helix nor strand, then skip
      if( !SS_TYPE->IsStructured())
      {
        return;
      }

      // this function will be called whenever an SSE is constructed from an aa sequence when the main geometry is
      // calculated, so this is probably not a terribly useful message
      BCL_MessageDbg
      (
        "SetToIdealConformation will strip off side chain atoms except CB for " + GetStaticClassName( *this)
      );

      storage::List< Atom>::iterator atom_itr( m_Atoms.Begin()), atom_itr_end( m_Atoms.End());

      // delete non backbone atoms
      while( atom_itr != atom_itr_end)
      {
        if( atom_itr->GetType()->IsBackBone() || atom_itr->GetType() == GetAtomTypes().CB)
        {
          // set the backbone to ideal conformation
          atom_itr->SetCoordinates
          (
            linal::Vector3D
            (
              atom_itr->GetType()->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
            ).Transform( TRANSFORMATION_MATRIX_3D)
          );

          // got to next atom
          ++atom_itr;
        }
        else
        {
          atom_itr = m_Atoms.Remove( atom_itr);
        }
      }

      // reinitialize the atom pointer vector
      m_AtomList = util::SiPtrVector< const Atom>( m_Atoms.Begin(), m_Atoms.End());
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator for AAComplete class
    //! @param AA_COMPLETE AAComplete object to be copied
    //! @return this after assignment to m_Atoms
    AAComplete &AAComplete::operator =( const AAComplete &AA_COMPLETE)
    {
      // assign base class AABase and data members
      AABase::operator =( AA_COMPLETE);

      // assign the atoms
      m_Atoms = AA_COMPLETE.m_Atoms;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AAComplete::Read( std::istream &ISTREAM)
    {
      // read base class
      AABase::Read( ISTREAM);

      // read member
      io::Serialize::Read( m_Atoms, ISTREAM);

      // reinitialize the siptrvector
      m_AtomList = util::SiPtrVector< const Atom>( m_Atoms.Begin(), m_Atoms.End());

      //return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &AAComplete::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write base class
      AABase::Write( OSTREAM, INDENT) << '\n';

      // write members
      io::Serialize::Write( m_Atoms, OSTREAM, INDENT);

      //return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief create a list of default backbone atoms, except CB (firstsidechainatom)
    //! @return List with  atoms, CA, C, N, O
    storage::List< Atom> AAComplete::DefaultBackboneAtoms()
    {
      // list of backbone atoms
      storage::List< Atom> backbone_atoms;

      // insert an atom for each backbone atom type
      for
      (
        storage::Set< AtomType>::const_iterator
          itr( GetAtomTypes().GetBackBoneAtomTypes().Begin()), itr_end( GetAtomTypes().GetBackBoneAtomTypes().End());
        itr != itr_end;
        ++itr
      )
      {
        backbone_atoms.PushBack( Atom( *itr));
      }

      // end
      return backbone_atoms;
    }

    //! @brief create a list of default backbone atoms
    //! @return List with  atoms, CA, C, N, O
    storage::List< Atom> AAComplete::DefaultBackboneAtomsWithFirstSideChainAtom( const AAType &AA_TYPE)
    {
      // get the backbone atoms
      storage::List< Atom> backbone_atoms( DefaultBackboneAtoms());

      // if type is Glycine
      if( AA_TYPE == GetAATypes().GLY)
      {
        backbone_atoms.PushBack( Atom( GetAtomTypes().HA2));
      }
      else
      {
        backbone_atoms.PushBack( Atom( GetAtomTypes().CB));
      }

      // end
      return backbone_atoms;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_aa.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_classes.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AA::s_Instance
    (
      GetObjectInstances().AddInstance( new AA())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AA::AA() :
      AABase()
    {
    }

    //! @brief construct AAData from util::ShPtr to AA
    //! @param SP_AA_DATA ShPtr to AAData to be copied
    AA::AA( const util::ShPtr< AAData> &SP_AA_DATA) :
      AABase( SP_AA_DATA)
    {
    }

    //! @brief constructor from AABase
    //! @param AA_BASE AABase
    AA::AA( const AABase &AA_BASE) :
      AABase( AA_BASE)
    {
    }

    //! @brief copy constructor - makes just a soft copy of m_Data
    //! @param AMINO_ACID AA to be copied
    AA::AA( const AA &AMINO_ACID) :
      AABase( AMINO_ACID)
    {
    }

    //! @brief virtual copy constructor
    AA *AA::Clone() const
    {
      return new AA( *this);
    }

    //! @brief virtual empty constructor with AAData
    //! @param SP_AA_DATA AAData object with information for the new AA
    //! this function is designed to be used in cases where AAClass is used for production multiple AA's
    //! which do not share one common AAData
    AA *AA::Empty( const util::ShPtr< AAData> &SP_AA_DATA) const
    {
      return new AA( SP_AA_DATA);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AA::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get types of atoms
    //! @return set of AtomTypes
    const storage::Set< AtomType> &AA::GetTypesOfAtoms() const
    {
      // initialize set of undefined atom types
      static const storage::Set< AtomType> s_atom_type_set;

      // return
      return s_atom_type_set;
    }

    //! @brief get all atoms
    //! @return SiPtrVector of Atoms
    const util::SiPtrVector< const Atom> &AA::GetAtoms() const
    {
      static const util::SiPtrVector< const Atom> s_undefined_atoms;
      return s_undefined_atoms;
    }

    //! @brief get the specified atom
    //! @brief ATOM_TYPE AtomType of interest
    //! @return atom with type ATOM_TYPE
    const Atom &AA::GetAtom( const AtomType &ATOM_TYPE) const
    {
      // initialize undefined atom
      static const Atom s_undefined_atom;

      // return
      return s_undefined_atom;
    }

    //! @brief get AAClass
    //! @return AAClass
    const AAClass &AA::GetAAClass() const
    {
      return GetAAClasses().e_AA;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator for AA class
    //! @param AA_RHS AA object to be copied
    //! @return this after assignment to AA_RHS is done
    AA &AA::operator =( const AA &AA_RHS)
    {
      // assign base class
      AABase::operator =( AA_RHS);

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AA::Read( std::istream &ISTREAM)
    {
      // read AABase
      AABase::Read( ISTREAM);

      //end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &AA::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write bases
      AABase::Write( OSTREAM, INDENT);

      //end
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_aa_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAData::s_Instance
    (
      GetObjectInstances().AddInstance( new AAData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct AAData from optional AAType and other information
    //! @param AA_TYPE amino acid type
    //! @param SEQ_ID sequence id
    //! @param PDB_ID pdb-file ID
    //! @param PDB_I_CODE insertion code for pdb-residues
    //! @param CHAIN_ID chain id of the amino acid
    AAData::AAData
    (
      const AAType &AA_TYPE,
      const int SEQ_ID,
      const int PDB_ID,
      const char PDB_I_CODE,
      const char CHAIN_ID
    ) :
      m_Type( AA_TYPE),
      m_SeqID( SEQ_ID),
      m_PdbID( PDB_ID),
      m_PdbICode( PDB_I_CODE),
      m_ChainID( CHAIN_ID),
      m_BlastProfile(),
      m_SSPredictions(),
      m_ExposurePrediction( util::GetUndefined< double>())
    {
    }

    //! @brief copy constructor
    AAData *AAData::Clone() const
    {
      return new AAData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return Type
    const AAType &AAData::GetType() const
    {
      return m_Type;
    }

    //! @brief return SeqID
    //! @return SeqID
    const int &AAData::GetSeqID() const
    {
      return m_SeqID;
    }

    //! @brief sets the sequence id to the given value
    //@ @param ID new sequence id
    void AAData::SetSeqID( int ID)
    {
      m_SeqID = ID;
    }

    //! @brief return PdbID
    //! @return PdbID
    const int &AAData::GetPdbID() const
    {
      return m_PdbID;
    }

    //! @brief return PdbICode
    //! @return PdbICode
    const char &AAData::GetPdbICode() const
    {
      return m_PdbICode;
    }

    //! @brief return chain id
    //! @return chain id
    char AAData::GetChainID() const
    {
      return m_ChainID;
    }

    //! @brief sets the chain id to the given value
    //! @param ID new chain id
    void AAData::SetChainID( char ID)
    {
      m_ChainID = ID;
    }

    //! @brief return blast profile
    //! @return blast profile
    const util::ShPtr< BlastProfile> &AAData::GetBlastProfile() const
    {
      return m_BlastProfile;
    }

    //! @brief return SSPredictions set
    //! @return SSPredictions set
    const storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> > &AAData::GetSSPredictions() const
    {
      return m_SSPredictions;
    }

    //! @brief return SSPredictions for given SS_METHOD
    //! @param SS_METHOD method of interest
    //! @return SSPredictions for given SS_METHOD
    util::SiPtr< const sspred::MethodInterface> AAData::GetSSPrediction( const sspred::Method &SS_METHOD) const
    {
      // search for the predictions for SS_METHOD in this map and store the iterator
      storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> >::const_iterator method_itr
      (
        m_SSPredictions.Find( SS_METHOD)
      );

      // if equal to end
      if( method_itr == m_SSPredictions.End())
      {
        // return undefined
        return util::SiPtr< const sspred::MethodInterface>();
      }

      // if found return itr dereferenced
      return method_itr->second;
    }

    //! @brief get the exposure prediction
    //! @return exposure prediction
    const double &AAData::GetExposurePrediction() const
    {
      return m_ExposurePrediction;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief smaller than, that compares the position in the sequence
    //! @param RHS_AA_DATA
    //! @return true if given AAData comes after this aa data by chain id, or for the same chain by seqid
    bool AAData::operator <( const AAData &RHS_AA_DATA) const
    {
      // compare chain id
      if( m_ChainID < RHS_AA_DATA.m_ChainID)
      {
        return true;
      }
      else if( m_ChainID > RHS_AA_DATA.m_ChainID)
      {
        return false;
      }

      // compare sequence id
      if( m_SeqID < RHS_AA_DATA.m_SeqID)
      {
        return true;
      }

      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AAData::Read( std::istream &ISTREAM)
    {
      // read data
      io::Serialize::Read( m_Type, ISTREAM);
      io::Serialize::Read( m_SeqID, ISTREAM);
      io::Serialize::Read( m_PdbID, ISTREAM);
      io::Serialize::Read( m_PdbICode, ISTREAM);
      io::Serialize::Read( m_ChainID, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &AAData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write data
      io::Serialize::Write( m_Type, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SeqID, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PdbID, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PdbICode, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ChainID, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_aa_sequence.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_compare.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_atom.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "util/bcl_util_sh_ptr_list.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AASequence::s_Instance
    (
      GetObjectInstances().AddInstance( new AASequence())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AASequence::AASequence() :
      m_Data(),
      m_ChainID( s_DefaultChainID),
      m_FastaHeader( GetDefaultFastaHeader())
    {
      Initialize( GetAAClasses().e_AA, 0);
    }

    //! @brief construct from AACLASS, LENGTH, CHAIN_ID and FASTA_HEADER
    //! @param AACLASS AAClass type for m_Data
    //! @param LENGTH number of residues in this sequence
    //! @param CHAIN_ID chain id of sequence
    //! @param FASTA_HEADER fasta header
    AASequence::AASequence
    (
      const AAClass &AACLASS,
      const size_t LENGTH,
      const char CHAIN_ID,
      const std::string &FASTA_HEADER
    ) :
      m_Data(),
      m_ChainID( CHAIN_ID)
    {
      SetFastaHeader( FASTA_HEADER);
      Initialize( AACLASS, LENGTH);
    }

    // TODO: can not ensure the people pass DATA from a sequence with AAClass of same type with this->m_AAClass
    //! @brief construct from ShPtrVector of AABases DATA, CHAIN_ID and FASTA_HEADER
    //! @param DATA ShPtrVector of AABases
    //! @param CHAIN_ID chain id of sequence
    //! @param FASTA_HEADER fasta header
    AASequence::AASequence
    (
      const util::ShPtrVector< AABase> &DATA,
      const char CHAIN_ID,
      const std::string &FASTA_HEADER
    ) :
      m_Data( DATA.HardCopy()),
      m_ChainID( CHAIN_ID),
      m_FastaHeader()
    {
      BCL_Assert( IsContinuous( m_Data), "given data are not continuous amino acids");
      BCL_Assert
      (
        m_Data.IsEmpty() || m_Data.FirstElement()->GetChainID() == m_ChainID,
        "first amino acid has different chain id than passed chainid! first amino acid chain id |" +
        m_Data.FirstElement()->GetIdentification() + "| The given chain id is  |" + m_ChainID + "|"
      );
      SetFastaHeader( FASTA_HEADER);
    }

    //! @brief copy constructor from another AA_SEQUENCE
    //! @param AA_SEQUENCE AASequence to be copied
    AASequence::AASequence( const AASequence &AA_SEQUENCE) :
      m_Data( AA_SEQUENCE.m_Data.HardCopy()),
      m_ChainID( AA_SEQUENCE.m_ChainID),
      m_FastaHeader( AA_SEQUENCE.m_FastaHeader)
    {
    }

    //! @brief virtual copy constructor
    AASequence *AASequence::Clone() const
    {
      return new AASequence( *this);
    }

    //! @brief hard copy all amino acids and their data
    //! @return aa sequence with independent hard copied AADatas
    AASequence *AASequence::HardCopy() const
    {
      // copy for this sequence
      AASequence *new_sequence( Clone());

      // iterate over all amino acids to make hard copy of the Data
      for
      (
        AASequence::iterator aa_itr( new_sequence->Begin()), aa_itr_end( new_sequence->End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // hard copy the data for that aa
        ( *aa_itr)->SetData( ( *aa_itr)->GetData().HardCopy());
      }

      // end
      return new_sequence;
    }

    //! @brief construct a new SequenceInterface from an AASequence interface implementation, sequence id and members
    //! @param ID sequence identifier
    //! @param MEMBERS sequence members
    //! @param CHAIN_ID chain id used to build the sequence
    //! @return pointer to a SequenceInterface
    align::SequenceInterface< AABase> *AASequence::Construct
    (
      const std::string &ID,
      const std::string &MEMBERS,
      const char CHAIN_ID
    ) const
    {
      AASequence *p_aa_sequence
      (
        AASequenceFactory::BuildSequenceFromFASTAString( MEMBERS, GetAAClasses().e_AA, CHAIN_ID).Clone()
      );
      p_aa_sequence->SetFastaHeader( ">" + ID);
      return p_aa_sequence;
    }

    //! @brief destructor
    AASequence::~AASequence()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AASequence::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Get default fasta header
    //! @return default fasta header
    const std::string &AASequence::GetDefaultFastaHeader()
    {
      // initialize static const default fasta header
      static const std::string s_default_fasta_header( "BCL_AASequence");

      // return
      return s_default_fasta_header;
    }

    //! @brief get identification of this AASequence
    //! @return string with identification for this AASequence
    std::string AASequence::GetSequenceIdentification() const
    {
      // initialize identification composed of chain ID, first seq id, the sequence one-letter codes, last seq-id
      // example: "chain: 'A' 23 VLALLV 28"
      std::string identification( "chain: " + util::Format()( m_ChainID));

      // add identifiers for the first and last residues this sse is spanning
      if( GetSize() > 0)
      {
        identification += ' ' + util::Format()( GetFirstAA()->GetSeqID())
          + ' ' + Sequence() + ' ' + util::Format()( GetLastAA()->GetSeqID());
      }

      // end
      return identification;
    }

    //! @brief add a AABase to the sequence with type based on ONE_LETTER_CODE
    //! @param ONE_LETTER_CODE char encoding the type of the AA
    void AASequence::AddMember( const char &ONE_LETTER_CODE)
    {
      // if m_Data is empty, set a default type of AAClass, otherwise get the type from the first AA
      const AAClass new_member_type( m_Data.IsEmpty() ? GetAAClasses().e_AA : GetFirstAA()->GetAAClass());

      // create ShPtr to AAData and set aatype, seq_id, chain_id
      util::ShPtr< AAData> new_member_data
      (
        new AAData
        (
          GetAATypes().AATypeFromOneLetterCode( toupper( ONE_LETTER_CODE)),
          m_Data.IsEmpty() ? AAData::s_DefaultSeqID : m_Data.LastElement()->GetSeqID() + 1,
          AAData::s_DefaultPdbID,
          AAData::s_DefaultPdbICode,
          m_ChainID
        )
      );

      // create ShPtr to AABase
      util::ShPtr< AABase> new_member( ( *new_member_type)->Empty( new_member_data));

      // push back in m_Data
      m_Data.PushBack( new_member);
    }

    //! @brief initializes the sequence with a new type of AAClass and fills with LENGTH number of residues
    //! deletes all existing aminoacids and builds up a new sequence with empty AABase interfaces
    //! @param AA_CLASS AAClass that identifies the amino acid type to be used
    //! @param LENGTH number of residues to be inserted
    void AASequence::Initialize( const AAClass &AA_CLASS, const size_t LENGTH)
    {
      // reset m_Data
      m_Data.Reset();

      // allocate memory for LENGTH number of residues
      m_Data.AllocateMemory( LENGTH);

      // iterate until LENGTH
      for( size_t i( 0); i < LENGTH; ++i)
      {
        util::ShPtr< AAData> sp_aa_data
        (
          new AAData( GetAATypes().e_Undefined, i + 1, i + 1, AAData::s_DefaultPdbICode, m_ChainID)
        );

        // create a ShPtr of AABase type with instance of AA_CLASS behind it
        util::ShPtr< AABase> sp_aa( ( *AA_CLASS)->Empty( sp_aa_data));

        // insert this new amino acid ShPtr into m_Data
        m_Data.PushBack( sp_aa);
      }
    }

    //! @brief get all atoms of all residues in this sequence
    //! @return SiPtrVector of Atoms of all residues in this sequence
    util::SiPtrVector< const Atom> AASequence::GetAtoms() const
    {
      // initialize atoms SiPtrVector to be returned
      util::SiPtrVector< const Atom> atoms;

      // allocate memory that is enough for total number of atoms to be inserted
      if( !m_Data.IsEmpty())
      {
        atoms.AllocateMemory( GetFirstAA()->GetNumberOfAtoms() * GetSize());
      }

      // iterate over all amino acids
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        //get atoms from the amino acid
        atoms.Append( ( *aa_itr)->GetAtoms());
      }

      // return
      return atoms;
    }

    //! @brief get all atoms for the specified atom types for all residues in the sequence
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return atoms of types specified in ATOM_TYPES  for all residues in the sequence
    util::SiPtrVector< const Atom> AASequence::GetAtoms
    (
      const storage::Set< AtomType> &ATOM_TYPES
    ) const
    {
      // initialize atoms SiPtrVector to be returned
      util::SiPtrVector< const Atom> atoms;

      // allocate memory that is enough for storing all atoms of specified ATOM_TYPES
      atoms.AllocateMemory( ATOM_TYPES.GetSize() * GetSize());

      // iterate over all amino acids
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // get atoms from the amino acid for specified ATOM_TYPES
        atoms.Append( ( *aa_itr)->GetAtoms( ATOM_TYPES));
      }

      // return
      return atoms;
    }

    //! @brief get all atom coordinates for all residues in the sequence
    //! @return all atom coordinates for all residues in the sequence
    util::SiPtrVector< const linal::Vector3D> AASequence::GetAtomCoordinates() const
    {
      // initialize coordinates SiPtrVector to be returned
      util::SiPtrVector< const linal::Vector3D> coordinates;

      // allocate memory that is enough for storing all atom coordinates
      if( !m_Data.IsEmpty())
      {
        coordinates.AllocateMemory( GetFirstAA()->GetNumberOfAtoms() * GetSize());
      }

      // iterate over all amino acids
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()),
         aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // get coordinates from the amino acid
        coordinates.Append( ( *aa_itr)->GetAtomCoordinates());
      }

      // return
      return coordinates;
    }

    //! @brief get all atom coordinates for the specified atom types for all residues in the sequence
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return coordinates of atoms of types specified in ATOM_TYPES  for all residues in the sequence
    util::SiPtrVector< const linal::Vector3D> AASequence::GetAtomCoordinates
    (
      const storage::Set< AtomType> &ATOM_TYPES
    ) const
    {
      // initialize coordinates SiPtrVector to be returned
      util::SiPtrVector< const linal::Vector3D> coordinates;

      // allocate memory that is enough for storing all atom coordinates of specified ATOM_TYPES
      coordinates.AllocateMemory( ATOM_TYPES.GetSize() * GetSize());

      // iterate over all amino acids
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // get coordinates from the amino acid for specified ATOM_TYPES
        coordinates.Append( ( *aa_itr)->GetAtomCoordinates( ATOM_TYPES));
      }

      // return
      return coordinates;
    }

    //! @brief find and return the iterator to the amino acid with the given sequence id
    //! @param SEQ_ID sequence id of the amino acid of interest
    //! @return the iterator to the amino acid with the given sequence id
    AASequence::iterator AASequence::FindAABySeqID( const int SEQ_ID)
    {
      // find and return by seq_id
      return std::find_if( Begin(), End(), AACompareBySeqID( SEQ_ID));
    }

    //! @brief find and return the iterator to the amino acid with the given sequence id
    //! @param SEQ_ID sequence id of the amino acid of interest
    //! @return the iterator to the amino acid with the given sequence id
    AASequence::const_iterator AASequence::FindAABySeqID( const int SEQ_ID) const
    {
      // find and return by seq_id
      return std::find_if( Begin(), End(), AACompareBySeqID( SEQ_ID));
    }

    //! @brief calculates the phi and psi for the amino acid with the given seq id
    //! @param SEQ_ID sequence id of the amino acid of interest
    //! @return pair of phi and psi values for the specified amino acid
    storage::VectorND< 2, double> AASequence::CalculatePhiPsi( const int SEQ_ID) const
    {
      // static undefined vector
      static const storage::VectorND< 2, double> s_undefined_vector( util::GetUndefined< double>());

      // find the amino acid
      const AASequence::const_iterator aa_itr( FindAABySeqID( SEQ_ID));

      // if not found return undefined
      if( aa_itr == m_Data.End())
      {
        return s_undefined_vector;
      }

      // construct the pair
      storage::VectorND< 2, double> phi_psi_pair( s_undefined_vector);

      // if the first residue
      if( aa_itr != m_Data.Begin())
      {
        phi_psi_pair.First() = ( *aa_itr)->CalculatePhi( ( *( aa_itr - 1))->GetAtom( GetAtomTypes().C));
      }
      else
      {
        phi_psi_pair.First() = ( *aa_itr)->Phi();
      }

      // if not the last residue
      if( aa_itr != m_Data.Last())
      {
        phi_psi_pair.Second() = ( *aa_itr)->CalculatePsi( ( *( aa_itr + 1))->GetAtom( GetAtomTypes().N));
      }
      else
      {
        phi_psi_pair.Second() = ( *aa_itr)->Psi();
      }

      // return
      return phi_psi_pair;
    }

    //! @brief set chain id of the sequence and all residues in the sequence to given CHAIN_ID
    //! @param CHAIN_ID new chain id of the sequence
    void AASequence::SetChainID( const char CHAIN_ID)
    {
      // first update the chain id of the sequence to CHAIN_ID
      m_ChainID = CHAIN_ID;

      // iterate over all residues in the sequence
      for
      (
        AASequence::iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // create a new amino acid with new aa data with the new CHAIN_ID
        util::ShPtr< AABase> sp_aa( ( *aa_itr)->Clone());
        util::ShPtr< AAData> sp_aa_data
        (
          new AAData( sp_aa->GetType(), sp_aa->GetSeqID(), sp_aa->GetPdbID(), sp_aa->GetPdbICode(), CHAIN_ID)
        );

        sp_aa->SetData( sp_aa_data);
        if( ( *aa_itr)->GetBlastProfilePtr().IsDefined())
        {
          sp_aa->SetBlastProfile( ( *aa_itr)->GetBlastProfile());
        }
        sp_aa->SetSSPredictions( ( *aa_itr)->GetSSPredictions());
        *aa_itr = sp_aa;
      }
    }

    //! @brief returns the amino acids with the given type
    //! @param AA_TYPE type of the amino acid to return
    //! @return shared pointers to the amino acids of the the given type
    const util::ShPtrList< AABase> AASequence::GetData( const AAType &AA_TYPE) const
    {
      // iterate over the sequence to find the amino acids of the given type
      util::ShPtrList< AABase> aa_list;
      typedef util::ShPtrVector< AABase>::const_iterator aa_it;
      for( aa_it it( Begin()), it_end( End()); it != it_end; ++it)
      {
        if( ( *it)->GetType() == AA_TYPE)
        {
          aa_list.Append( *it);
        }
      }

      return aa_list;
    }

    //! @brief change the fasta header of the sequence to given FASTA_HEADER
    //! @param FASTA_HEADER new fasta header of the sequence
    void AASequence::SetFastaHeader( const std::string &FASTA_HEADER)
    {
      m_FastaHeader = util::TrimString( FASTA_HEADER);

      if( !m_FastaHeader.empty() && m_FastaHeader[ 0] == s_FastaHeaderChar)
      {
        m_FastaHeader = m_FastaHeader.substr( 1, m_FastaHeader.size() - 1);
      }
      m_FastaHeader = util::TrimString( m_FastaHeader);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION translation vector to be applied
    void AASequence::Translate( const linal::Vector3D &TRANSLATION)
    {
      //iterate over all amino acids
      for
      (
        AASequence::iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // translate the amino acid coordinates
        ( *aa_itr)->Translate( TRANSLATION);
      }
    }

    //! @brief transform the object by a given TRANSFORMATION_MATRIX_3D
    //! @param TRANSFORMATIONMATRIX3D TransformationMatrix3D to be applied
    void AASequence::Transform( const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D)
    {
      //iterate over all amino acids
      for
      (
        AASequence::iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // transform the amino acid coordinates
        ( *aa_itr)->Transform( TRANSFORMATIONMATRIX3D);
      }
    }

    //! @brief rotate the object by a given ROTATION_MATRIX_3D
    //! @param ROTATIONMATRIX3D RotationMatrix3D to be applied
    void AASequence::Rotate( const math::RotationMatrix3D &ROTATIONMATRIX3D)
    {
      //iterate over all amino acids
      for
      (
        AASequence::iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // rotate the amino acid coordinates
        ( *aa_itr)->Rotate( ROTATIONMATRIX3D);
      }
    }

    //! return util::ShPtrVector of AASequences chopped to pieces of size SIZE with one aa gap
    util::ShPtrVector< AASequence> AASequence::ChopSequence( const size_t SIZE) const
    {
      BCL_Assert( util::IsDefined( SIZE), "undefined size supplied");

      //calculate the numer of pieces
      size_t number_of_sequences = ( GetSize() + 1) / ( SIZE + 1);

      //instantiate util::ShPtrVector to pieces
      util::ShPtrVector< AASequence> chopped_sequence;

      //nothing to be chopped
      if( number_of_sequences < 2)
      {
        chopped_sequence.PushBack( util::ShPtr< AASequence>( new AASequence( *this)));
        return chopped_sequence;
      }

      // pushback last sequence which does not provide a gap
      chopped_sequence.PushBack
                       (
                         util::ShPtr< AASequence>
                         (
                           new AASequence
                           (
                             SubSequence
                             (
                               0, //position of subsequence
                               size_t( math::Absolute( double( 1) * double( GetSize()) / double( number_of_sequences))) //length of subsequence
                             )
                           )
                         )
                       );

      //collect chopped sequences and make an one aagap between sequences
      for( size_t i( 1); i < number_of_sequences - 1; ++i)
      {
        chopped_sequence.PushBack
        (
          util::ShPtr< AASequence>
          (
            new AASequence
            (
              SubSequence
              (
                size_t( math::Absolute( double( i    ) * double( GetSize()) / double( number_of_sequences))) + 1,  //position of subsequence
                size_t( math::Absolute( double( i + 1) * double( GetSize()) / double( number_of_sequences)))       //length of subsequence
                - size_t( math::Absolute( double( i    ) * double( GetSize()) / double( number_of_sequences))) - 1 //length of subsequence
              )
            )
          )
        );
      }

      // pushback last sequence which does not provide a gap
      chopped_sequence.PushBack
      (
        util::ShPtr< AASequence>
        (
          new AASequence
          (
            SubSequence
            (
              size_t( math::Absolute( double( number_of_sequences - 1) * double( GetSize()) / double( number_of_sequences))) + 1, //position of subsequence
              GetSize() - size_t( math::Absolute( double( number_of_sequences - 1) * double( GetSize()) / double( number_of_sequences))) - 1 //length of subsequence
            )
          )
        )
      );

      //return chopped sequences
      return chopped_sequence;
    }

    //! @brief clips number of amino acidss from the beginning and the end of the sequence
    //! @param NUMBER_AMINO_ACIDS number of amino acids to be clipped from each end
    //! @return this sequence after clipping NUMBER_AMINO_ACIDS from each end
    AASequence &AASequence::ClipEnds( const size_t NUMBER_AMINO_ACIDS)
    {
      // check that after clipping we are left with at least 1 amino acid in the sequence
      if( m_Data.GetSize() <= 2 * NUMBER_AMINO_ACIDS + 1)
      {
        BCL_MessageCrt
        (
          "sequence too short: " + util::Format()( m_Data.GetSize()) +
          " <= " + util::Format()( 2 * NUMBER_AMINO_ACIDS + 1)
        );

        // end
        return *this;
      };
      // create a subsequence accordingly
      operator =( SubSequence( NUMBER_AMINO_ACIDS, m_Data.GetSize() - 2 * NUMBER_AMINO_ACIDS));

      // return
      return *this;
    }

    //! @brief prepends new amino acid to the beginning of the sequence
    //! @param AMINO_ACID amino acid to be prepended to the sequence
    void AASequence::PushFront( const AABase &AMINO_ACID)
    {
      if( m_Data.IsEmpty())
      {
        m_ChainID = AMINO_ACID.GetChainID();
      }
      else
      {
        // check that amino acid precedes sequence
        BCL_Assert
        (
          AMINO_ACID.DoesPrecede( *GetFirstAA()),
          "amino acid does not preced this: " + AMINO_ACID.GetIdentification() + " .. " + GetSequenceIdentification()
        );
      }

      // clone the AMINO_ACID and pushback it into the m_Data
      m_Data.InsertElement( m_Data.Begin(), util::ShPtr< AABase>( AMINO_ACID.Clone()));
    }

    //! @brief appends util::ShPtr to amino acid to the end of the sequence
    //! @param SP_AMINO_ACID ShPtr to amino acid to be appended to the end of the sequence
    void AASequence::PushBack( const util::ShPtr< AABase> &SP_AMINO_ACID)
    {
      if( m_Data.IsEmpty())
      {
        m_ChainID = SP_AMINO_ACID->GetChainID();
      }
      else
      {
        // check that amino acids follows this sequence
        BCL_Assert
        (
          GetLastAA()->DoesPrecede( *SP_AMINO_ACID),
          "amino acid does not follow this: " + GetSequenceIdentification() + " .. " + SP_AMINO_ACID->GetIdentification()
        );
      }

      // push back SP_AMINO_ACID
      m_Data.PushBack( SP_AMINO_ACID);
    }

    //! @brief check if a second sequence follows this by seqid and chain id (this precedes given sequence)
    //! @param SEQUENCE AASequence that should follow this
    //! @return true if neither seq is empty and this last aa precedes arguments first aa
    bool AASequence::DoesPrecede( const AASequence &SEQUENCE) const
    {
      return
           !m_Data.IsEmpty() && !SEQUENCE.m_Data.IsEmpty() // at least one aa in each data
        && m_Data.LastElement()->DoesPrecede( *SEQUENCE.m_Data.FirstElement()); // checks seqid and chainid
    }

    //! @brief appends SEQUENCE to the end of this sequence
    //! @param SEQUENCE AASequence to be appended to this sequence
    void AASequence::AppendSequence( const AASequence &SEQUENCE)
    {
      if( SEQUENCE.m_Data.IsEmpty())
      {
        return;
      }

      if( m_Data.IsEmpty())
      {
        m_ChainID = SEQUENCE.GetChainID();
        m_Data = SEQUENCE.m_Data.HardCopy();
      }
      else
      {
        // assure that chain ids do match
        BCL_Assert
        (
          DoesPrecede( SEQUENCE),
          "given sequence does not follow this: " + GetSequenceIdentification() +
          " .. " + SEQUENCE.GetSequenceIdentification()
        );
        // insert the aminoacids of the SEQUENCE
        m_Data.InsertElements( m_Data.GetSize(), SEQUENCE.m_Data.HardCopy());
      }
    }

    //! @brief prepends SEQUENCE to the end of this sequence
    //! @param SEQUENCE AASequence to be prepended to this sequence
    void AASequence::PrependSequence( const AASequence &SEQUENCE)
    {
      if( SEQUENCE.m_Data.IsEmpty())
      {
        return;
      }

      if( m_Data.IsEmpty())
      {
        m_ChainID = SEQUENCE.GetChainID();
        m_Data = SEQUENCE.m_Data.HardCopy();
      }
      else
      {
        // assure that chain ids do match
        BCL_Assert
        (
          SEQUENCE.DoesPrecede( *this),
          "given sequence does precede this: " + SEQUENCE.GetSequenceIdentification() +
          " .. " + GetSequenceIdentification()
        );

        // insert the amino acids of the SEQUENCE
        m_Data.InsertElements( 0, SEQUENCE.m_Data.HardCopy());
      }
    }

    //! @brief return sequence as string
    //! @return sequence as a string
    std::string AASequence::Sequence() const
    {
      // instantiate string for one letter code sequence
      std::string sequence;

      // iterate over all amino acids in the sequence
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // concatenate one letter code for this amino acid
        sequence += ( *aa_itr)->GetType()->GetOneLetterCode();
      }

      // end
      return sequence;
    }

    //! @brief let the aadata of the given sequence point to the aadata of this sequence - determined by pdbID
    //! @param AA_SEQUENCE sequence of amino acids which aas will be connected to the aadata of this sequence
    //! @param SET_ATOMS set to true if atoms should also be connected by this function
    void AASequence::ConnectAADataByPdbID
    (
      AASequence &AA_SEQUENCE,
      const bool &SET_ATOMS
    ) const
    {
      // if given sequence is empty, return
      if( AA_SEQUENCE.m_Data.IsEmpty())
      {
        // set the sequence's chain id to this chain id
        AA_SEQUENCE.m_ChainID = m_ChainID;

        return;
      }

      // instantiate iterator to the begin and end of this sequence and given AA_SEQUENCE
      AASequence::const_iterator aa_itr( m_Data.Begin()), aa_itr_end( m_Data.End());
      AASequence::iterator arg_aa_itr( AA_SEQUENCE.Begin()), arg_aa_itr_end( AA_SEQUENCE.End());
      BCL_Assert( !m_Data.IsEmpty(), "m_Data.IsEmpty()");
      BCL_Assert( !AA_SEQUENCE.m_Data.IsEmpty(), "AA_SEQUENCE.IsEmpty()");

      // find first matching amino acid
      while
      (
        // iterate as long PdbID or ICode are not equal
        ( *aa_itr)->GetPdbID() != ( *arg_aa_itr)->GetPdbID() ||
        ( *aa_itr)->GetPdbICode() != ( *arg_aa_itr)->GetPdbICode()
      )
      {
        // move further in this sequence until the match is found
        ++aa_itr;

        // if the end has been reached
        if( aa_itr == aa_itr_end)
        {
          // issue warning and return
          BCL_MessageCrt( "could not find any matching starting amino acid");
          return;
        }
      }

      // iterate over the residues in the portion of the full sequence that is equal to the given sequence
      while( aa_itr != aa_itr_end && arg_aa_itr != arg_aa_itr_end)
      {
        // update the data of the residue in AASEQUENCE to corresponding residues data in this sequence
        ( *arg_aa_itr)->SetData( ( *aa_itr)->GetData());
        if( SET_ATOMS)
        {
          ( **arg_aa_itr).SetAtoms( ( **aa_itr).GetAtoms());
        }

        // move further in both sequences
        ++aa_itr;
        ++arg_aa_itr;

        // if end has been reached at either of the sequence
        if( aa_itr == aa_itr_end && arg_aa_itr != arg_aa_itr_end)
        {
          // exit
          BCL_Exit( "there are not enough amino acids in this sequence to be connected to the aas in the argument", -1);
        }
      }

      // set the sequence's chain id to this chain id
      AA_SEQUENCE.m_ChainID = m_ChainID;

      // end
      return;
    }

    //! @brief returns whether all coordinates for the given atom types are defined
    //! @param ATOM_TYPES Atom Types of interest
    bool AASequence::HasDefinedCoordinates
    (
      const storage::Set< AtomType> &ATOM_TYPES// = biol::GetAtomTypes().GetBackBoneAtomTypes()
    ) const
    {
      // iterate over the residues
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()), aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // if there are undefined coords
        if( !( *aa_itr)->HasDefinedCoordinates())
        {
          return false;
        }
      }

      // this point is reached only if all coordinates were defined, therefore return true
      return true;
    }

    //! @brief counts the number of residues that have defined coordinates for all the given atom types
    //! @param ATOM_TYPES Atom Types of interest
    size_t AASequence::CountDefinedAACoordinates
    (
      const storage::Set< AtomType> &ATOM_TYPES // = GetAtomTypes().GetBackBoneAtomTypes()
    ) const
    {
      size_t number_defined( 0);
      // iterate over the residues
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()), aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // if there are undefined coords
        if( ( *aa_itr)->HasDefinedCoordinates( ATOM_TYPES))
        {
          ++number_defined;
        }
      }

      // this point is reached only if all coordinates were defined, therefore return true
      return number_defined;
    }

  ///////////////
  // operators //
  ///////////////

    //! operator = for assigning this sequence to AA_SEQUENCE
    //! @param AA_SEQUENCE AASequence to be copied
    //! @return this sequence after being assigned to AA_SEQUENCE
    AASequence &AASequence::operator =( const AASequence &AA_SEQUENCE)
    {
      // assign data members
      m_Data = AA_SEQUENCE.m_Data.HardCopy();
      m_ChainID = AA_SEQUENCE.m_ChainID;
      m_FastaHeader = AA_SEQUENCE.m_FastaHeader;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write fasta to given OSTREAM in blocks of BLOCK_SIZE
    //! @param OSTREAM output stream
    //! @param BLOCK_SIZE number of amino acids to be outputted in each block
    //! @return std::ostream which was written to
    std::ostream &AASequence::WriteFasta( std::ostream &OSTREAM, const size_t BLOCK_SIZE) const
    {
      // output the fasta header
      OSTREAM << s_FastaHeaderChar << m_FastaHeader << '\n';

      // initialize counter
      size_t counter( 0);

      // iterate over the amino acids in the sequence
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()), aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // output the one letter code of this amino acid
        ( *aa_itr)->WriteFasta( OSTREAM);

        // if BLOCK_SIZE has been reached
        if( !( ( counter + 1) % BLOCK_SIZE))
        {
          // output endline
          OSTREAM << '\n';
        }

        // increment counter
        ++counter;
      }

      // line break at end of output
      OSTREAM << '\n';

      // end
      return OSTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AASequence::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);
      io::Serialize::Read( m_ChainID, ISTREAM);
      io::Serialize::Read( m_FastaHeader, ISTREAM);
      SetFastaHeader( m_FastaHeader);

      // check that chain id is correct and the amino acids are continuous
      BCL_Assert
      (
           m_Data.IsEmpty()
        || ( IsContinuous( m_Data) && ( m_Data.FirstElement()->GetChainID() == m_ChainID)),
        "read sequence is not continuous"
      );

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &AASequence::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ChainID    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FastaHeader, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief check if a given vector of amino acids is continuous
    //! @param DATA amino acids
    //! @return true if two consecutive amino acids follow each other by chain id and seq id
    bool AASequence::IsContinuous( const util::ShPtrVector< AABase> &DATA)
    {
      if( DATA.IsEmpty())
      {
        return true;
      }

      // iterate over amino acids
      const_iterator itr( DATA.Begin()), itr_follow( itr), itr_end( DATA.End());
      ++itr_follow;
      for( ; itr_follow != itr_end; ++itr, ++itr_follow)
      {
        if( !( *itr)->DoesPrecede( **itr_follow))
        {
          return false;
        }
      }

      return true;
    }

    //! checks whether two sequences have overlap
    bool DoOverlap( const AASequence &SEQUENCE_A, const AASequence &SEQUENCE_B)
    {
      //check for identical chain IDs
      if( SEQUENCE_A.GetChainID() != SEQUENCE_B.GetChainID())
      {
        BCL_MessageDbg( "different chains");
        return false;
      }

      // empty sequences do not overlaps
      if( SEQUENCE_A.GetSize() == 0 || SEQUENCE_B.GetSize() == 0)
      {
        return false;
      }

      //check if the range between beginning and ending seqid of each sequence are not overlapping
      return SEQUENCE_A.GetLastAA()->GetSeqID() >= SEQUENCE_B.GetFirstAA()->GetSeqID()
             && SEQUENCE_A.GetFirstAA()->GetSeqID() <= SEQUENCE_B.GetLastAA()->GetSeqID();
    }

    //! @brief calculates the distance in sequence between SEQUENCE_A and SEQUENCE_B
    //! @param SEQUENCE_A first sequence
    //! @param SEQUENCE_B second sequence
    //! @return the distance between SEQUENCE_A and SEQUENCE_B, Undefined if they overlap or are in different chains
    size_t CalculateSequenceDistance( const AASequence &SEQUENCE_A, const AASequence &SEQUENCE_B)
    {
      // if different chains or sequences overlap return undefined size_t
      if( SEQUENCE_A.GetChainID() != SEQUENCE_B.GetChainID() || DoOverlap( SEQUENCE_A, SEQUENCE_B))
      {
        return util::GetUndefinedSize_t();
      }
      // if SEQUENCE_A comes first
      if( SEQUENCE_A.GetLastAA()->GetSeqID() < SEQUENCE_B.GetFirstAA()->GetSeqID())
      {
        // return the distance from last residue of SEQUENCE_A to first residue of SEQUENCE_B
        return SEQUENCE_B.GetFirstAA()->GetSeqID() - SEQUENCE_A.GetLastAA()->GetSeqID() - 1;
      }
      // else if SEQUENCE_B comes first
      else
      {
        // return the distance from last residue of SEQUENCE_A to first residue of SEQUENCE_B
        return SEQUENCE_A.GetFirstAA()->GetSeqID() - SEQUENCE_B.GetLastAA()->GetSeqID() - 1;
      }
    }

    //! @brief gives the distance between the c atom of the n terminal aa and the n atom of the c terminal aa
    //! @param N_TERMINAL_AA the aa that provides the c in the peptide bond
    //! @param C_TERMINAL_AA the aa that provides the n in the peptide bond
    //! @return distance between the c atom of the n terminal aa and the n atom of the c terminal aa
    double GetPeptideBondLength( const AABase &N_TERMINAL_AA, const AABase &C_TERMINAL_AA)
    {
      if
      (
        !N_TERMINAL_AA.GetAtom( GetAtomTypes().C).GetCoordinates().IsDefined() ||
        !C_TERMINAL_AA.GetAtom( GetAtomTypes().N).GetCoordinates().IsDefined()
      )
      {
        return util::GetUndefinedDouble();
      }

      return linal::Distance
        (
          N_TERMINAL_AA.GetAtom( GetAtomTypes().C).GetCoordinates(),
          C_TERMINAL_AA.GetAtom( GetAtomTypes().N).GetCoordinates()
        );
    }

    //! @brief generates subsequences of specified radius around every residue and returns them
    //! @param SEQUENCE AASequence of interest
    //! @param WINDOW_RADIUS radius of the window
    //! @param UNDEFINED_AA undefined amino acid to be used for windows of border amino acids
    //! @return map of window sequences
    storage::Map< util::SiPtr< const AABase>, util::SiPtrVector< const AABase> > CreateWindowsFromAminoAcids
    (
      const AASequence &SEQUENCE,
      const size_t WINDOW_RADIUS,
      const AABase &UNDEFINED_AA
    )
    {
      // initialize window map
      storage::Map< util::SiPtr< const AABase>, util::SiPtrVector< const AABase> > window_map;

      // initialize SiPtr to undefined AA
      util::SiPtr< const AABase> undefined_aa( UNDEFINED_AA);

      // initialize window length
      const size_t window_length( 2 * WINDOW_RADIUS + 1);

      // make sure sequence is long enough for at least one window
      if( SEQUENCE.GetSize() <= window_length)
      {
        BCL_MessageStd
        (
          "The given sequence's length is smaller than the window radius " + util::Format()( SEQUENCE.GetSize()) +
            " < " + util::Format()( WINDOW_RADIUS)
         );
        return window_map;
      }

      // initialize amino acid counter
      size_t aa_count( 0);

      // now iterate over residues
      for
      (
        AASequence::const_iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr, ++aa_count
      )
      {
        // initialize the SiPtrVector and store the reference for this amino acid
        util::SiPtrVector< const AABase> &current_window( window_map[ **aa_itr]);

        // if there are not enough residues on the left of this residue fill with undefined residues
        if( aa_count < WINDOW_RADIUS)
        {
          // calculate how many undefined residues are different
          const size_t number_undefined( WINDOW_RADIUS - aa_count);

          // insert number_undefined counts of given UNDEFINED_AA
          current_window.Append
          (
            storage::Vector< util::SiPtr< const AABase> >( number_undefined, undefined_aa)
          );

          // append the rest of the the sequence
          current_window.Append
          (
            SEQUENCE.GetData().SubShPtrVector( 0, ( window_length - number_undefined))
          );
        }
        // if there are not enough residues on the right side of this residue fill with undefined residues
        else if( SEQUENCE.GetSize() <= WINDOW_RADIUS + aa_count)
        {
          // calculate how many undefined residues are different
          const size_t number_undefined( aa_count + WINDOW_RADIUS + 1 - SEQUENCE.GetSize());

          // append the rest of the the sequence
          current_window.Append
          (
            SEQUENCE.GetData().SubShPtrVector( aa_count - WINDOW_RADIUS, ( window_length - number_undefined))
          );

          // insert number_undefined counts of given UNDEFINED_AA
          current_window.Append
          (
            storage::Vector< util::SiPtr< const AABase> >( number_undefined, undefined_aa)
          );
        }
        // else it is safe to make a simple subsequence
        else
        {
          // append the subsequence
          current_window.Append
          (
            SEQUENCE.GetData().SubShPtrVector( aa_count - WINDOW_RADIUS, window_length)
          );
        }
      }

      // end
      return window_map;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_aa_sequence_factory.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_aa_back_bone_completer.h"
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "biol/bcl_biol_aa_sequence_phi_psi.h"
#include "fold/bcl_fold_mutate_aa_set_phi.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically
#include <iterator>

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AASequenceFactory::s_Instance
    (
      GetObjectInstances().AddInstance( new AASequenceFactory())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AASequenceFactory::AASequenceFactory()
    {
    }

    //! @brief Clone function
    //! @return pointer to new AASequenceFactory
    AASequenceFactory *AASequenceFactory::Clone() const
    {
      return new AASequenceFactory( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AASequenceFactory::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief appends a residue to the C-terminal end of a sequence
    //! @param SEQUENCE AASequence to be extended
    //! @param AMINO_ACID to be appended
    //! @param PHI phi angle for added residue
    //! @param PSI psi angle for last residue in original sequence
    void AASequenceFactory::AppendAA
    (
      AASequence &SEQUENCE,
      const AABase &AMINO_ACID,
      const double PHI,
      const double PSI
    )
    {
      AppendSequence
      (
        SEQUENCE,
        AASequence
        (
          util::ShPtrVector< AABase>( 1, util::ShPtr< AABase>( AMINO_ACID.Clone())),
          SEQUENCE.GetChainID(),
          SEQUENCE.GetFastaHeader()
        ),
        PHI,
        PSI
      );
    }

    //! @brief appends a sequence to the C-terminal end of a sequence
    //! @param N_TERMINAL_SEQUENCE AASequence to be extended
    //! @param C_TERMINAL_SEQUENCE AASequence to be appended
    //! @param PHI phi angle for first residue in C-terminal sequence
    //! @param PSI psi angle for last residue in N-terminal sequence
    void AASequenceFactory::AppendSequence
    (
      AASequence &N_TERMINAL_SEQUENCE,
      const AASequence &C_TERMINAL_SEQUENCE,
      const double PHI,
      const double PSI
    )
    {
      // change psi angle
      {
        const int seq_id( N_TERMINAL_SEQUENCE.GetLastAA()->GetSeqID());
        const storage::VectorND< 2, double> new_phi_psi( 0.0, PSI);
        storage::VectorND< 2, double> new_phi_psi_change( AASequenceFlexibility::CalculatePhiPsiChange( N_TERMINAL_SEQUENCE, seq_id, new_phi_psi));
        new_phi_psi_change.First() = 0;
        AASequenceFlexibility::ChangePhiPsi( N_TERMINAL_SEQUENCE, seq_id, new_phi_psi_change, AASequenceFlexibility::e_CTerminal);
      }

      // calculate necessary transformation, apply and prepend
      const math::TransformationMatrix3D trans_c_term( TransformationAppend( N_TERMINAL_SEQUENCE, *C_TERMINAL_SEQUENCE.GetFirstAA(), PHI));
      AASequence c_term_seq( C_TERMINAL_SEQUENCE);
      c_term_seq.Transform( trans_c_term);
      N_TERMINAL_SEQUENCE.AppendSequence( c_term_seq);
    }

    //! @brief transformation to append a c term sequence onto the given n-term sequence
    //! @param N_TERMINAL_SEQUENCE AASequence fixed in space
    //! @param C_TERMINAL_AA amino acid to be appended
    //! @param PHI angle for first residue in cterm
    math::TransformationMatrix3D AASequenceFactory::TransformationAppend
    (
      const AASequence &N_TERMINAL_SEQUENCE,
      const AABase &C_TERMINAL_AA,
      const double PHI
    )
    {
      // make sure that the cterminal sequence has its coordinates defined
      BCL_Assert
      (
        N_TERMINAL_SEQUENCE.HasDefinedCoordinates(),
        "N-terminal aa sequence does not have all defined coordinates " + util::Format()( N_TERMINAL_SEQUENCE)
      );

      // make sure the n-terminal sequence has its coordinates defined
      BCL_Assert
      (
        C_TERMINAL_AA.HasDefinedCoordinates(),
        "C-terminal aa does not have all defined coordinates " + util::Format()( C_TERMINAL_AA)
      );

      // get the coordinates of the c atom of the last residue in the n-terminal sequence
      const linal::Vector3D &anchor_aa_coord_c
      (
        N_TERMINAL_SEQUENCE.GetLastAA()->GetAtom( GetAtomTypes().C).GetCoordinates()
      );

      // overall transformation
      math::TransformationMatrix3D overall_trans;

      // create ShPtr to a new c-terminal sequence cloned from "C_TERMINAL_SEQUENCE"
      util::ShPtr< AABase> new_c_terminal_aa( C_TERMINAL_AA.Clone());

    ///////////////////
    // superimpose N //
    ///////////////////

      {
        // calculate the coordinates of where the N of first residue in the cterminal sequence should go according to
        // "PSI"
        const Atom desired_nitrogen
        (
          AABackBoneCompleter::GenerateN( *N_TERMINAL_SEQUENCE.GetLastAA())
        );

        // get the current position of the nitrogen of the first residue of the c-terminal sequence
        const linal::Vector3D &current_coord_n
        (
          new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates()
        );

        // translate c-terminal sequence to correct position (at desired N position)
        const linal::Vector3D c_super_trans( desired_nitrogen.GetCoordinates() - current_coord_n);
        new_c_terminal_aa->Translate( c_super_trans);

        // accumulate overall transformation
        overall_trans( c_super_trans);
      }

    ///////////////////
    // superimpose C //
    ///////////////////

      {
        // desired n
        const Atom desired_carbon( AABackBoneCompleter::GenerateC( *new_c_terminal_aa, PHI));

        // cross product
        const linal::Vector3D rotation_axis
        (
          linal::CrossProduct
          (
            new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates(),
            desired_carbon.GetCoordinates(),
            anchor_aa_coord_c
          )
        );

        // rotation angle
        const double angle
        (
          linal::ProjAngle
          (
            new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates(),
            desired_carbon.GetCoordinates(),
            anchor_aa_coord_c
          )
        );

        // transformation matrix
        math::TransformationMatrix3D transformation( -new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates());

        // apply rotation to superimpose n onto n
        transformation( math::RotationMatrix3D( rotation_axis, -angle));

        // move back
        transformation( new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates());

        // apply to amino acid and overall transformation
        new_c_terminal_aa->Transform( transformation);
        overall_trans( transformation);
      }

    ///////////////
    // set omega //
    ///////////////

      {
        // calculate the rotation necessary in order to properly set omega to 180 degrees
        const double omega_rotation
        (
          AABackBoneCompleter::s_OmegaCACNCA - // optimal omega
          new_c_terminal_aa->CalculateOmega( N_TERMINAL_SEQUENCE.GetLastAA()->GetCA(), N_TERMINAL_SEQUENCE.GetLastAA()->GetAtom( GetAtomTypes().C))
        );

        // create transformation matrix and add the necessary transformations to set the omega angle
        math::TransformationMatrix3D omega_transform
        (
          -new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates()
        );

        // add the rotation to "c_n_ca_transform" to set the omega angle
        omega_transform
        (
          math::RotationMatrix3D
          (
            new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates() - anchor_aa_coord_c,
            -omega_rotation
          )
        );

        // add the translation to move the c-terminal sequence back to its original position but rotated
        omega_transform( new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates());

        // now transform the c-terminal sequence so that omega is properly set
        new_c_terminal_aa->Transform( omega_transform);

        // accumulate overall transformation
        overall_trans( omega_transform);
      }

    /////////////
    // set phi //
    /////////////

      {
        // transform the c-terminal sequence so that the phi of the first residue is set correctly to "PHI"
        const math::TransformationMatrix3D phi_trans
        (
          fold::MutateAASetPhi
          (
            N_TERMINAL_SEQUENCE.GetLastAA()->GetAtom( GetAtomTypes().C), PHI
          ).GetTransformationMatrix( *new_c_terminal_aa)
        );
//        new_c_terminal_aa->Transform
//        (
//          phi_trans
//        );

        // accumulate overall transformation
        overall_trans( phi_trans);
      }

      // return the overall transformation
      return overall_trans;
    }

    //! @brief prepends a residue to the N-terminal end of a sequence
    //! @param AMINO_ACID to be prepended
    //! @param SEQUENCE AASequence to be extended
    //! @param PHI phi angle for first residue in original sequence
    //! @param PSI psi angle for added residue
    void AASequenceFactory::PrependAA
    (
      const AABase &AMINO_ACID,
      AASequence &SEQUENCE,
      const double PHI,
      const double PSI
    )
    {
      PrependSequence
      (
        AASequence
        (
          util::ShPtrVector< AABase>( 1, util::ShPtr< AABase>( AMINO_ACID.Clone())),
          SEQUENCE.GetChainID(),
          SEQUENCE.GetFastaHeader()
        ),
        SEQUENCE,
        PHI,
        PSI
      );
    }

    //! @brief prepends a sequence to the N-terminal end of a sequence
    //! @param N_TERMINAL_SEQUENCE AASequence to be prepended
    //! @param C_TERMINAL_SEQUENCE AASequence to be extended
    //! @param PHI phi angle for first residue in C-terminal sequence
    //! @param PSI psi angle for last residue in N-terminal sequence
    void AASequenceFactory::PrependSequence
    (
      const AASequence &N_TERMINAL_SEQUENCE,
      AASequence &C_TERMINAL_SEQUENCE,
      const double PHI,
      const double PSI
    )
    {
      // create ShPtr to a new n-terminal sequence cloned from "N_TERMINAL_SEQUENCE"
      AASequence new_n_terminal_sequence( N_TERMINAL_SEQUENCE);

      // change psi angle
      {
        const int seq_id( new_n_terminal_sequence.GetLastAA()->GetSeqID());
        const storage::VectorND< 2, double> new_phi_psi( 0.0, PSI);
        storage::VectorND< 2, double> new_phi_psi_change( AASequenceFlexibility::CalculatePhiPsiChange( new_n_terminal_sequence, seq_id, new_phi_psi));
        new_phi_psi_change.First() = 0;
        AASequenceFlexibility::ChangePhiPsi( new_n_terminal_sequence, seq_id, new_phi_psi_change, AASequenceFlexibility::e_NTerminal);
      }

      // transformation matrix
      const math::TransformationMatrix3D trans_n_term( TransformationPrepend( *new_n_terminal_sequence.GetLastAA(), C_TERMINAL_SEQUENCE, PHI));
      new_n_terminal_sequence.Transform( trans_n_term);
      C_TERMINAL_SEQUENCE.PrependSequence( new_n_terminal_sequence);
    }

    //! @brief transformation to append a c term sequence onto the given n-term sequence
    //! @param N_TERMINAL_AA amino acid to be appended
    //! @param C_TERMINAL_SEQUENCE AASequence fixed in space
    //! @param PHI angle for first residue in cterm
    math::TransformationMatrix3D AASequenceFactory::TransformationPrepend
    (
      const AABase &N_TERMINAL_AA,
      const AASequence &C_TERMINAL_SEQUENCE,
      const double PHI
    )
    {
      // make sure that the n-terminal residue has its coordinates defined
      BCL_Assert
      (
        N_TERMINAL_AA.HasDefinedCoordinates(),
        "N-terminal aa does not have all defined coordinates " + util::Format()( N_TERMINAL_AA)
      );

      // make sure that the c-terminal sequence has its coordinates defined
      BCL_Assert
      (
        C_TERMINAL_SEQUENCE.HasDefinedCoordinates(),
        "C-terminal aa sequence does not have all defined coordinates " + util::Format()( C_TERMINAL_SEQUENCE)
      );

      // get the coordinates of the n atom of the first residue in the c-terminal sequence
      const linal::Vector3D &anchor_aa_coord_n
      (
        C_TERMINAL_SEQUENCE.GetFirstAA()->GetAtom( GetAtomTypes().N).GetCoordinates()
      );

      // create ShPtr to a new n-terminal aa
      util::ShPtr< AABase> new_n_term_aa( N_TERMINAL_AA.Clone());

      // overall transformation
      math::TransformationMatrix3D overall_trans;

    ///////////////////
    // superimpose C //
    ///////////////////

      {
        // calculate the coordinates of where the C of last residue in the nterminal sequence should go according to
        // "PHI"
        const Atom desired_carbon
        (
          AABackBoneCompleter::GenerateC( *C_TERMINAL_SEQUENCE.GetFirstAA(), PHI)
        );

        // get the current position of the carbon of the last residue of the n-terminal sequence
        const linal::Vector3D &current_coord_c
        (
          new_n_term_aa->GetAtom( GetAtomTypes().C).GetCoordinates()
        );

        // translate n-terminal sequence to correct position (at desired C position)
        const linal::Vector3D super_c_trans( desired_carbon.GetCoordinates() - current_coord_c);
        new_n_term_aa->Translate( super_c_trans);

        // accumulate overall transformation
        overall_trans( super_c_trans);
      }

    ///////////////////
    // superimpose N //
    ///////////////////

      {
        // desired n
        const Atom desired_nitrogen( AABackBoneCompleter::GenerateN( *new_n_term_aa));

        // cross product
        const linal::Vector3D rotation_axis
        (
          linal::CrossProduct
          (
            new_n_term_aa->GetAtom( GetAtomTypes().C).GetCoordinates(),
            desired_nitrogen.GetCoordinates(),
            anchor_aa_coord_n
          )
        );

        // rotation angle
        const double angle
        (
          linal::ProjAngle
          (
            new_n_term_aa->GetAtom( GetAtomTypes().C).GetCoordinates(),
            desired_nitrogen.GetCoordinates(),
            anchor_aa_coord_n
          )
        );

        // transformation matrix
        math::TransformationMatrix3D transformation( -new_n_term_aa->GetAtom( GetAtomTypes().C).GetCoordinates());

        // apply rotation to superimpose n onto n
        transformation( math::RotationMatrix3D( rotation_axis, -angle));

        // move back
        transformation( new_n_term_aa->GetAtom( GetAtomTypes().C).GetCoordinates());

        // apply to amino acid and overall transformation
        new_n_term_aa->Transform( transformation);
        overall_trans( transformation);
      }

    ///////////////
    // set omega //
    ///////////////

      {
        // calculate the rotation necessary in order to properly set omega to 180 degrees
        const double omega_rotation
        (
          AABackBoneCompleter::s_OmegaCACNCA - // optimal omega
          C_TERMINAL_SEQUENCE.GetFirstAA()->CalculateOmega( new_n_term_aa->GetCA(), new_n_term_aa->GetAtom( GetAtomTypes().C))
        );

        // create transformation matrix and add the necessary transformations to set the omega angle
        math::TransformationMatrix3D omega_transform
        (
          -anchor_aa_coord_n
        );

        // creat rotation matrix around peptide bond  and add to omega tranformation
        omega_transform
        (
          math::RotationMatrix3D
          (
            anchor_aa_coord_n - new_n_term_aa->GetAtom( GetAtomTypes().C).GetCoordinates(),
            omega_rotation //< rotate in the correct direction
          )
        );

        // add the translation to move the n-terminal sequence back to its original position but rotated
        omega_transform( anchor_aa_coord_n);

        // now transform the n-terminal sequence so that omega is properly set
        new_n_term_aa->Transform( omega_transform);

        // accumulate overall transformation
        overall_trans( omega_transform);
      }

    ////////////////////////
    // phi transformation //
    ////////////////////////

      {
        // transform the n-terminal aa so that the phi is correct relative to the c terminal sequence
        const math::TransformationMatrix3D phi_trans
        (
          fold::MutateAASetPhi
          (
            new_n_term_aa->GetAtom( GetAtomTypes().C), PHI
          ).GetTransformationMatrix( *C_TERMINAL_SEQUENCE.GetFirstAA())
        );
//        new_n_term_aa->Transform
//        (
//          phi_trans
//        );

        // accumulate overall transformation
        overall_trans( phi_trans);
      }

      // return the overall transformation applied
      return overall_trans;
    }

    //! @brief fits the passed sequence to the phi/psi angles
    //! @param SEQUENCE AASequence to be used
    //! @param PHI_PSI phi and psi values to be applied
    //! @param SS_TYPE sstype to be applied if the sequence is longer than the given phi/psi information
    void AASequenceFactory::FitSequence
    (
      AASequence &SEQUENCE,
      const AASequencePhiPsi &PHI_PSI,
      const SSType &SS_TYPE
    )
    {
      // if the sequence is empty or no phi/psi's are given
      if( SEQUENCE.GetSize() == 0 || PHI_PSI.GetAngles().IsEmpty())
      {
        // return with no changes
        return;
      }

      // idealize the sequence according to the given type
      IdealizeSequence( SEQUENCE, SS_TYPE);

      // get the middle index
      const size_t middle_seq_index( SEQUENCE.GetSize() / 2);

      // move the sequence so that the middle residues are superimposed
      SEQUENCE.Transform
      (
        quality::RMSD::SuperimposeCoordinates
        (
          util::SiPtrVector< const linal::Vector3D>::Create
          (
            PHI_PSI.GetN(),
            PHI_PSI.GetCA(),
            PHI_PSI.GetC()
          ),
          util::SiPtrVector< const linal::Vector3D>::Create
          (
            SEQUENCE.GetData()( middle_seq_index)->GetAtom( GetAtomTypes().N).GetCoordinates(),
            SEQUENCE.GetData()( middle_seq_index)->GetAtom( GetAtomTypes().CA).GetCoordinates(),
            SEQUENCE.GetData()( middle_seq_index)->GetAtom( GetAtomTypes().C).GetCoordinates()
          )
        )
      );

      // determine the corresponding indeces to be used
      const size_t middle_angle_index( PHI_PSI.GetAngles().GetSize() / 2);
      size_t first_index( 0);
      size_t last_index( PHI_PSI.GetAngles().GetSize() - 1);
      int first_seq_id( SEQUENCE.GetFirstAA()->GetSeqID());

      // if the sequence has more residues than phi/psi angles given
      if( SEQUENCE.GetSize() > PHI_PSI.GetAngles().GetSize())
      {
        // move up the first seq id to match with the beginning of the angle vector
        first_seq_id += ( SEQUENCE.GetSize() - PHI_PSI.GetAngles().GetSize()) / 2;
      }

      // if the angle vector is larger than the sequence
      if( PHI_PSI.GetAngles().GetSize() > SEQUENCE.GetSize())
      {
        // adjust the indeces to match the vector size
        first_index += ( PHI_PSI.GetAngles().GetSize() - SEQUENCE.GetSize()) / 2;
        last_index = first_index + SEQUENCE.GetSize() - 1;

        // move the first seq id to compensate
        first_seq_id -= first_index;
      }

      // iterate backwards from the residue
      for( int i( middle_angle_index - 1); i >= int( first_index); --i)
      {
        // set the previous phi/psi's
        AASequenceFlexibility::SetPhiPsi
        (
          SEQUENCE, first_seq_id + i, PHI_PSI.GetAngles()( i), AASequenceFlexibility::e_NTerminal
        );
      }

      // set the middle residue
      AASequenceFlexibility::SetPhiPsi
      (
        SEQUENCE,
        first_seq_id + middle_angle_index,
        PHI_PSI.GetAngles()( middle_angle_index),
        AASequenceFlexibility::e_Bidirectional
      );

      // iterate forwards from the residue
      for( int i( middle_angle_index + 1); i <= int( last_index); ++i)
      {
        // set the next phi/psi's
        AASequenceFlexibility::SetPhiPsi
        (
          SEQUENCE, first_seq_id + i, PHI_PSI.GetAngles()( i), AASequenceFlexibility::e_CTerminal
        );
      }
    }

    //! @brief Idealize the coordinates for the given AASequence according to given SSType
    //! @param AA_SEQUENCE AASequence to be idealized
    //! @param SS_TYPE SSType to be used for idealization
    void AASequenceFactory::IdealizeSequence( AASequence &AA_SEQUENCE, const SSType &SS_TYPE)
    {
      // find the transformation matrix that transforms one residue to the following residue in this SSType,
      // including the translation related to z axis( half the height of this SSE)
      math::TransformationMatrix3D transform
      (
        linal::Vector3D( 0.0, 0.0, -1.0 * AA_SEQUENCE.GetSize() * 0.5 * SS_TYPE->GetRiseInZPerResidue())
      );

      // iterate over aa's
      for
      (
        AASequence::iterator
          aa_itr( AA_SEQUENCE.GetData().Begin()), aa_itr_end( AA_SEQUENCE.GetData().End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // first set the residues to ideal transformed conformation
        ( *aa_itr)->SetToIdealConformation( SS_TYPE, transform);

        // transform the transformation matrix
        transform( SS_TYPE->GetTransformationMatrixForResidues());
      }
    }

    //! @brief read fasta from given ISTREAM and initialize the sequence with instances of AA_CLASS
    //! @param ISTREAM input stream
    //! @param AA_CLASS type of amino acid to be used
    //! @param CHAIN_ID chain id to be used for the sequence
    //! @return newly constructed AASequence
    AASequence AASequenceFactory::BuildSequenceFromFASTA
    (
      std::istream &ISTREAM,
      const AAClass &AA_CLASS,
      const char CHAIN_ID
    )
    {
      // local variables
      char c;

      // read header
      // whitespace is skipped automatically (unless you add "<< noskipws")
      std::string fasta_header( ">");
      ISTREAM >> c;

      // if ">" is found as first character, assume there is a header
      if( c == AASequence::s_FastaHeaderChar)
      {
        // read the complete first line and save in m_FastaHeader
        std::string line;
        std::getline( ISTREAM, line);
        fasta_header += line;
      }
      // if there is no ">"
      else
      {
        // put back the read character and create a new header
        ISTREAM.putback( c);
        fasta_header += AASequence::GetDefaultFastaHeader();
      }

      // initialize iterators on the stream to read the sequence
      std::istream_iterator< char> stream_itr( ISTREAM);
      std::istream_iterator< char> stream_itr_end;
      size_t id = 1;

      // initialize vector of aabase
      util::ShPtrVector< AABase> aa_data;

      // while there are chars left in STREAM
      while( stream_itr != stream_itr_end)
      {
        // if current_oneletter_code is "*", we are at the end of the sequence, then break
        if( *stream_itr == '*')
        {
          break;
        }

        // if current_oneletter_code is ">", we found the next sequence, so putback the char and end this sequence
        if( *stream_itr == AASequence::s_FastaHeaderChar)
        {
          ISTREAM.unget(); // unget the ">" character that the next sequence can be read correctly
          break;
        }

        // get aatype from current_oneletter_code
        const AAType current_aatype( GetAATypes().AATypeFromOneLetterCode( toupper( *stream_itr)));

        // if current aa type is undefined
        if
        (
          current_aatype == GetAATypes().e_Undefined &&
          ( *stream_itr) != GetAATypes().e_Undefined->GetOneLetterCode()
        )
        {
          // give user a message on every unkown aatype we found
          BCL_MessageCrt
          (
            std::string( "Unknown one letter code: ") + std::string( 1, *stream_itr)
            + " (ASCII: " + util::Format()( int( *stream_itr)) + ") supplied!"
          );
          break;
        }

        util::ShPtr< AAData> sp_aa_data( new AAData( current_aatype, int( id), int( id), AAData::s_DefaultPdbICode, CHAIN_ID));

        // initialize ShPtr to an amino acid of type AA_CLASS corresponding AAData
        util::ShPtr< AABase> sp_aa( ( *AA_CLASS)->Empty( sp_aa_data));

        // insert this aminoacid into the sequence
        aa_data.PushBack( sp_aa);

        // move further in the stream and in the sequence ids
        ++stream_itr;
        ++id;
      }

      // end
      return AASequence( aa_data, CHAIN_ID, fasta_header);
    }

    //! @brief read fasta from provided SEQUENCE_STRING and initialize the sequences with amino acids of type AA_CLASS
    //! @param SEQUENCE_STRING sequence string that has the one letter code sequence
    //! @param AA_CLASS type of amino acid to be used
    //! @param CHAIN_ID chain id to be used for the sequence
    //! @return AA sequence built from the string
    AASequence AASequenceFactory::BuildSequenceFromFASTAString
    (
      const std::string &SEQUENCE_STRING,
      const AAClass &AA_CLASS,
      const char CHAIN_ID
    )
    {
      // read the string into a stream
      std::stringstream ss;
      ss << SEQUENCE_STRING;

      // call BuildSequenceFromFASTA function with this stream
      return BuildSequenceFromFASTA( ss, AA_CLASS, CHAIN_ID);
    }

    //! @brief read sequence data, fasta and pdb files, from the given filenames and returns their sequences
    //! @param FILENAMES a vector of sequence filenames
    //! @return a ShPtrVector of sequences
    util::ShPtrVector< AASequence>
    AASequenceFactory::ReadSequenceData( const storage::Vector< std::string> &FILENAMES)
    {
      util::ShPtrVector< AASequence> sequences; // instantiate ShPtrVector to collect all sequences
      io::IFStream read; // IFStream for reading files
      pdb::Factory factory; // pdb::Factory to reading pdbs

      for
      (
        storage::Vector< std::string>::const_iterator itr( FILENAMES.Begin()), itr_end( FILENAMES.End());
        itr != itr_end;
        ++itr
      )
      {
        io::File::MustOpenIFStream( read, *itr);
        std::string extension( io::File::GetFullExtension( io::File::RemovePath( *itr)));

        // read sequences from different file formats
        // searching for the format extension in all extensions allow also compressed files
        if( extension.find( "fasta") != std::string::npos)
        {
          while( !read.eof()) // read all sequences in this file
          {
            // new AASequence to read next sequence
            util::ShPtr< AASequence> sequence( new AASequence( AASequenceFactory::BuildSequenceFromFASTA( read)));
            sequences.PushBack( sequence); // save sequence
          }
        }
        else if( extension.find( pdb::GetDefaultFileExtension()) != std::string::npos)
        {
          pdb::Handler pdb( read);
          util::ShPtrVector< AASequence> pdb_sequences( factory.AASequencesFromPDB( pdb));
          sequences.InsertElements( sequences.End(), pdb_sequences);
        }

        io::File::CloseClearFStream( read);
      }

      return sequences;
    }

    //! @brief calculates the transformation needed for superimposing SEQUENCE_A to SEQUENCE_B
    //! @param AA_SEQUENCE_A AASequence which will be superimposed to AA_SEQUENCE_B
    //! @param AA_SEQUENCE_B AASequence to which superimposition will be calculated against
    //! @return transformation matrix that gives the best superimposition
    math::TransformationMatrix3D AASequenceFactory::CalculateSuperimposition
    (
      const AASequence &AA_SEQUENCE_A,
      const AASequence &AA_SEQUENCE_B
    )
    {
      // if this is the same AASequence then return identity matrix
      if( &AA_SEQUENCE_A == &AA_SEQUENCE_B)
      {
        return math::TransformationMatrix3D();
      }

      // calculate transformation matrix for super imposition
      math::TransformationMatrix3D trans_matrix
      (
        quality::RMSD::SuperimposeCoordinates
        (
          AA_SEQUENCE_B.GetAtomCoordinates( GetAtomTypes().GetBackBoneAtomTypes()),
          AA_SEQUENCE_A.GetAtomCoordinates( GetAtomTypes().GetBackBoneAtomTypes())
        )
      );

      // problem with superimposition? try to superimpose all coordinates
      if( trans_matrix == math::TransformationMatrix3D())
      {
        trans_matrix = quality::RMSD::SuperimposeCoordinates
        (
          AA_SEQUENCE_B.GetAtomCoordinates(),
          AA_SEQUENCE_A.GetAtomCoordinates()
        );
      }

      // end
      return trans_matrix;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AASequenceFactory::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AASequenceFactory::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_aa_sequence_flexibility.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_atom.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! @brief conversion to a string from a SequenceDirection
    //! @param SEQUENCE_DIRECTION the sequence direction to get a string for
    //! @return a string representing that sequence direction
    const std::string &AASequenceFlexibility::GetSequenceDirectionName
    (
      const AASequenceFlexibility::SequenceDirection &SEQUENCE_DIRECTION
    )
    {
      static const std::string s_descriptors[] =
      {
        "n_terminal",
        "c_terminal",
        "bidirectional",
        GetStaticClassName< SequenceDirection>()
      };
      return s_descriptors[ size_t( SEQUENCE_DIRECTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new AASequenceFlexibility
    AASequenceFlexibility *AASequenceFlexibility::Clone() const
    {
      return new AASequenceFlexibility( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AASequenceFlexibility::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief convenience function to calculate the number of different phi/psi between two SSEs
    //! @param SEQUENCE_A first sequence of interest
    //! @param SEQUENCE_B second sequence  of interest
    //! @return number of different phi/psi pairs between given sequences
    size_t AASequenceFlexibility::GetNumberDifferentPhiPsi
    (
      const AASequence &SEQUENCE_A,
      const AASequence &SEQUENCE_B
    )
    {
      // make sure they are the same SSEs
      BCL_Assert
      (
        SEQUENCE_A.GetSize() == SEQUENCE_B.GetSize(),
        "The sequences have to be same size for comparison!\nSEQUENCE_A " +
          util::Format()( SEQUENCE_A.GetSize()) + "\t" + SEQUENCE_A.Sequence() + "\nvs\nSEQUENCE_B " +
          util::Format()( SEQUENCE_B.GetSize())
      );

      // initialize count
      size_t count( 0);

      const int seqid_begin( SEQUENCE_A.GetFirstAA()->GetSeqID());
      const int seqid_end( SEQUENCE_B.GetLastAA()->GetSeqID());

      // iterate over seqids
      for( int seqid( seqid_begin); seqid <= seqid_end; ++seqid)
      {
        // calculate phi psi from both
        const storage::VectorND< 2, double> phi_psi_a( SEQUENCE_A.CalculatePhiPsi( seqid));
        const storage::VectorND< 2, double> phi_psi_b( SEQUENCE_B.CalculatePhiPsi( seqid));

        // if first residue do not check
        const bool phi_match
        (
          seqid == seqid_begin ? true : math::EqualWithinTolerance( phi_psi_a.First(), phi_psi_b.First())
        );

        const bool psi_match
        (
          seqid == seqid_end ? true : math::EqualWithinTolerance( phi_psi_a.Second(), phi_psi_b.Second())
        );

        // increment if any failed
        count += !( phi_match && psi_match);
      }

      // end
      return count;
    }

    //! @brief function to calculate a Transformation from begin atom, end atom and the rotation angle
    //! @param ATOM_BEGIN Atom that forms the beginning point of the rotation axis
    //! @param ATOM_END Atom that forms the end point of the rotation axis
    //! @param ROTATION degrees of rotation
    math::TransformationMatrix3D AASequenceFlexibility::CalculateTransformation
    (
      const linal::Vector3D &ATOM_BEGIN,
      const linal::Vector3D &ATOM_END,
      const double ROTATION
    )
    {
      // initialize the transformation so that you move the ATOM_BEGIN to origin
      math::TransformationMatrix3D transform( -ATOM_BEGIN);

      // rotate around the vector connecting ATOM_BEGIN to ATOM_END
      transform( math::RotationMatrix3D( linal::Vector3D( ATOM_END - ATOM_BEGIN), ROTATION));

      // move ATOM_BEGIN to its original position
      transform( ATOM_BEGIN);

      // end
      return transform;
    }

    //! @brief convenience function to rotate a specific atom of a residue with given transformation
    //! @param AMINO_ACID Amino acid of interest
    //! @param ATOM_TYPE AtomType of interest
    //! @param TRANSFORMATION Transformation to be applied
    void AASequenceFlexibility::TransformAtom
    (
      AABase &AMINO_ACID,
      const AtomType &ATOM_TYPE,
      const math::TransformationMatrix3D &TRANSFORMATION
    )
    {
      // get the correponding atom from this amino acid
      const Atom &atom( AMINO_ACID.GetAtom( ATOM_TYPE));

      // if the type is defined, thus the atom of the type ATOM_TYPE exists
      if( atom.GetType().IsDefined())
      {
        // make a copy of the atom and apply the psi transformation and put it back
        Atom new_atom( atom);
        new_atom.Transform( TRANSFORMATION);
        AMINO_ACID.SetAtom( new_atom);
      }
    }

    //! @brief calculate the difference in phi and psi to be applied
    //! @param SEQUENCE Sequence to be bent
    //! @param SEQ_ID sequence id of the amino acid of interest
    //! @param PHI_PSI phi and psi angles to be set
    //! @return change in phi and psi angles to be applied
    storage::VectorND< 2, double> AASequenceFlexibility::CalculatePhiPsiChange
    (
      AASequence &SEQUENCE,
      const int SEQ_ID,
      const storage::VectorND< 2, double> &PHI_PSI
    )
    {
      // calculate the current phi, psi values for the given amino acid
      const storage::VectorND< 2, double> current_phi_psi( SEQUENCE.CalculatePhiPsi( SEQ_ID));

      // calculate the differences
      return
        storage::VectorND< 2, double>
        (
          PHI_PSI.First() - current_phi_psi.First(), PHI_PSI.Second() - current_phi_psi.Second()
        );
    }

    //! @brief set the phi psi of the given amino acid to given values
    //! @param SEQUENCE Sequence to be bent
    //! @param SEQ_ID sequence id of the amino acid of interest
    //! @param PHI_PSI phi and psi angles to be set
    //! @param SEQUENCE_DIRECTION enumerator whether the bending should be applied towards N-terminal, C-terminal or both directions
    //! @return transformations applied towards N-terminal and C-Terminal
    storage::VectorND< 2, math::TransformationMatrix3D> AASequenceFlexibility::SetPhiPsi
    (
      AASequence &SEQUENCE,
      const int SEQ_ID,
      const storage::VectorND< 2, double> &PHI_PSI,
      const AASequenceFlexibility::SequenceDirection &SEQUENCE_DIRECTION
    )
    {
      // calculate the change and call the corresponding change phi psi function
      return ChangePhiPsi( SEQUENCE, SEQ_ID, CalculatePhiPsiChange( SEQUENCE, SEQ_ID, PHI_PSI), SEQUENCE_DIRECTION);
    }

    //! @brief change the phi psi of the given amino acid by the given values
    //! @param SEQUENCE Sequence to be bent
    //! @param SEQ_ID sequence id of the amino acid of interest
    //! @param PHI_PSI_CHANGE change in phi and psi angles to be applied
    //! @param SEQUENCE_DIRECTION enumerator whether the bending should be applied towards N-terminal, C-terminal or both directions
    //! @return transformations applied towards N-terminal and C-Terminal
    storage::VectorND< 2, math::TransformationMatrix3D> AASequenceFlexibility::ChangePhiPsi
    (
      AASequence &SEQUENCE,
      const int SEQ_ID,
      const storage::VectorND< 2, double> &PHI_PSI_CHANGE,
      const AASequenceFlexibility::SequenceDirection &SEQUENCE_DIRECTION
    )
    {
      // get the phi and psi change
      double phi_change( PHI_PSI_CHANGE.First());
      const double psi_change( PHI_PSI_CHANGE.Second());

      // if PHI or PSI is undefined return
      if( !util::IsDefined( phi_change) && !util::IsDefined( psi_change))
      {
        BCL_MessageDbg
        (
          "The provided PHI:" + util::Format()( phi_change) +
          " and PSI: " + util::Format()( psi_change) + " is undefined!"
        );
        return
          storage::VectorND< 2, math::TransformationMatrix3D>( util::GetUndefined< math::TransformationMatrix3D>());
      }

      // find the residue
      const AASequence::iterator aa_itr( SEQUENCE.FindAABySeqID( SEQ_ID));

      // make sure the AA is found
      BCL_Assert
      (
        aa_itr != SEQUENCE.End(),
        "No Amino acid with seqid " + util::Format()( SEQ_ID) + " in sequence " + SEQUENCE.GetSequenceIdentification()
      );

      // proline can not be rotated around phi
      if( ( *aa_itr)->GetType() == GetAATypes().PRO)
      {
        phi_change = 0.0;
      }

      // store reference to coordinates of nitrogen, ca and c for this amino acid
      const linal::Vector3D &coord_N( ( *aa_itr)->GetAtom( GetAtomTypes().N).GetCoordinates());
      const linal::Vector3D &coord_CA( ( *aa_itr)->GetCA().GetCoordinates());
      const linal::Vector3D &coord_C( ( *aa_itr)->GetAtom( GetAtomTypes().C).GetCoordinates());

      // initialize the transformations for phi and psi
      math::TransformationMatrix3D phi_transform, psi_transform;

      // if phi change is non-zero
      if( util::IsDefined( phi_change) && phi_change != double( 0.0))
      {
        // if towards N terminal or both directions
        if( SEQUENCE_DIRECTION == e_NTerminal || SEQUENCE_DIRECTION == e_Bidirectional)
        {
          // calculate phi transformation, the change has to be negated since we rotate residues beforehand
          phi_transform = CalculateTransformation( coord_N, coord_CA, phi_change);
          TransformAtom( **aa_itr, GetAtomTypes().H, phi_transform);
        }
        // else towards C terminal
        else
        {
          // calculate the phi transformation
          phi_transform = CalculateTransformation( coord_N, coord_CA, -phi_change);

          // transform this amino acid, keep hydrogen in place
          const Atom hydrogen( ( *aa_itr)->GetAtom( GetAtomTypes().H));
          ( *aa_itr)->Transform( phi_transform);
          if( hydrogen.GetType().IsDefined())
          {
            ( *aa_itr)->SetAtom( hydrogen);
          }
        }
      }

      // if psi change is non-zero
      if( util::IsDefined( psi_change) && psi_change != double( 0.0))
      {
        // if towards C terminal or both directions
        if( SEQUENCE_DIRECTION == e_CTerminal || SEQUENCE_DIRECTION == e_Bidirectional)
        {
          // calculate psi transformation and transform oxygen for this residue
          psi_transform = CalculateTransformation( coord_CA, coord_C, -psi_change);
          TransformAtom( **aa_itr, GetAtomTypes().O, psi_transform);
        }
        // else towards N terminal
        else
        {
          // calculate psi transformation the angle has to be negated since we rotate residues towards N-terminal
          psi_transform = CalculateTransformation( coord_CA, coord_C, psi_change);
          // transform the residue but keep oxygen in place
          const Atom oxygen( ( *aa_itr)->GetAtom( GetAtomTypes().O));
          ( *aa_itr)->Transform( psi_transform);
          // reset oxygen
          ( *aa_itr)->SetAtom( oxygen);
        }
      }

      // calculate the total transformation from phi and psi
      math::TransformationMatrix3D total_transformation( phi_transform);
      total_transformation( psi_transform);

      // construct the return transformation applied to N-terminal and C-terminal
      storage::VectorND< 2, math::TransformationMatrix3D> transformation_pair;

      // switch over bending types
      switch( SEQUENCE_DIRECTION)
      {
        case e_NTerminal:
        {
          // update the transformation pair
          transformation_pair.First() = total_transformation;

          // iterate over the N terminal part of the sequence
          for( AASequence::iterator itr( SEQUENCE.Begin()); itr != aa_itr; ++itr)
          {
            // apply the total transformation to each amino acid
            ( *itr)->Transform( total_transformation);
          }
          break;
        }
        case e_CTerminal:
        {
          // update the transformation pair
          transformation_pair.Second() = total_transformation;

          // iterate over the C terminal part of the sequence
          for( AASequence::iterator itr( aa_itr + 1), itr_end( SEQUENCE.End()); itr != itr_end; ++itr)
          {
            // apply the total transformation to each amino acid
            ( *itr)->Transform( total_transformation);
          }
          break;
        }
        case e_Bidirectional:
        {
          // update the transformation pair
          transformation_pair.First() = phi_transform;
          transformation_pair.Second() = psi_transform;

          // iterate over the N terminal part of the sequence
          for( AASequence::iterator itr( SEQUENCE.Begin()); itr != aa_itr; ++itr)
          {
            // apply the total transformation to each amino acid
            ( *itr)->Transform( phi_transform);
          }

          // iterate over the C terminal part of the sequence
          for( AASequence::iterator itr( aa_itr + 1), itr_end( SEQUENCE.End()); itr != itr_end; ++itr)
          {
            // apply the total transformation to each amino acid
            ( *itr)->Transform( psi_transform);
          }
          break;
        }
        case s_NumberSequenceDirections:
        {
          break;
        }
      }

      // return the pair of transformations applied
      return transformation_pair;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AASequenceFlexibility::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AASequenceFlexibility::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_aa_sequence_phi_psi.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AASequencePhiPsi::s_Instance
    (
      GetObjectInstances().AddInstance( new AASequencePhiPsi())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AASequencePhiPsi::AASequencePhiPsi() :
      m_N(),
      m_CA(),
      m_C(),
      m_Angles()
    {
    }

    //! @brief construct from an AASequence
    //! @param SEQUENCE sequence to be used
    AASequencePhiPsi::AASequencePhiPsi( const AASequence &SEQUENCE) :
      m_N(),
      m_CA(),
      m_C(),
      m_Angles()
    {
      InitializeAngles( SEQUENCE);
    }

    //! @brief construct from the 3 coordinates and the phi/psi angles
    //! @param N_COORDS coordinates for N atom of middle residue
    //! @param CA_COORDS coordinates for CA atom of middle residue
    //! @param C_COORDS coordinates for C atom of middle residue
    //! @param ANGLES vector of phi-psi angles
    AASequencePhiPsi::AASequencePhiPsi
    (
      const linal::Vector3D &N_COORDS,
      const linal::Vector3D &CA_COORDS,
      const linal::Vector3D &C_COORDS,
      const storage::Vector< storage::VectorND< 2, double> > &ANGLES
    ) :
      m_N( N_COORDS),
      m_CA( CA_COORDS),
      m_C( C_COORDS),
      m_Angles( ANGLES)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AASequencePhiPsi
    AASequencePhiPsi *AASequencePhiPsi::Clone() const
    {
      return new AASequencePhiPsi( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AASequencePhiPsi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the geometric center of the object
    //! @return the geometric center of the object
    linal::Vector3D AASequencePhiPsi::GetCenter() const
    {
      return m_CA;
    }

    //! @brief returns whether the object is defined (first phi and last psi can be nan but everything else cannot)
    //! @return whether the object is defined
    bool AASequencePhiPsi::IsDefined() const
    {
      // if any of the coords are undefined
      if( !m_N.IsDefined() || !m_CA.IsDefined() || !m_C.IsDefined())
      {
        // return false
        return false;
      }

      // check if the first element psi value is undefined
      if( !util::IsDefined( m_Angles.FirstElement().Second()))
      {
        return false;
      }

      // check if the last element phi value is undefined
      if( !util::IsDefined( m_Angles.LastElement().First()))
      {
        return false;
      }

      // if the # of angles is greater than 2
      if( m_Angles.GetSize() > 2)
      {
        // get iterators
        storage::Vector< storage::VectorND< 2, double> >::const_iterator angle_itr( m_Angles.Begin());
        storage::Vector< storage::VectorND< 2, double> >::const_iterator angle_itr_end( m_Angles.End());

        // already checked first and last elements so skip them
        ++angle_itr;
        --angle_itr_end;

        // iterate over the remaining angles
        for( ; angle_itr != angle_itr_end; ++angle_itr)
        {
          // if either value is undefined
          if( !util::IsDefined( angle_itr->First()) || !util::IsDefined( angle_itr->Second()))
          {
            // return false
            return false;
          }
        }
      }

      // if this point is reached all values are defined
      return true;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION Translation to be applied
    void AASequencePhiPsi::Translate( const linal::Vector3D &TRANSLATION)
    {
      m_N.Translate( TRANSLATION);
      m_CA.Translate( TRANSLATION);
      m_C.Translate( TRANSLATION);
    }

    //! @brief transform the object by a given TransformationMatrix3D
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void AASequencePhiPsi::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      m_N.Transform( TRANSFORMATION_MATRIX_3D);
      m_CA.Transform( TRANSFORMATION_MATRIX_3D);
      m_C.Transform( TRANSFORMATION_MATRIX_3D);
    }

    //! @brief rotate the object by a given RotationMatrix3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
    void AASequencePhiPsi::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      m_N.Rotate( ROTATION_MATRIX_3D);
      m_CA.Rotate( ROTATION_MATRIX_3D);
      m_C.Rotate( ROTATION_MATRIX_3D);
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AASequencePhiPsi::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_N, ISTREAM);
      io::Serialize::Read( m_CA, ISTREAM);
      io::Serialize::Read( m_C, ISTREAM);
      io::Serialize::Read( m_Angles, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AASequencePhiPsi::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_N, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_C, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Angles, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief initializes member variables from an AASequence
    //! @param SEQUENCE sequence to be used
    void AASequencePhiPsi::InitializeAngles( const AASequence &SEQUENCE)
    {
      // if the sequence is empty or not continuous
      if( SEQUENCE.GetData().IsEmpty())
      {
        // warn the user and return
        BCL_MessageStd
        (
          "Given sequence is empty, so no phi-psi angles will be determined: " +
            SEQUENCE.GetSequenceIdentification()
        );
        return;
      }

      // get the middle index
      const size_t middle_index( SEQUENCE.GetSize() / 2);

      // get the coordinates of the first amino acid
      m_N = SEQUENCE.GetData()( middle_index)->GetAtom( GetAtomTypes().N).GetCoordinates();
      m_CA = SEQUENCE.GetData()( middle_index)->GetAtom( GetAtomTypes().CA).GetCoordinates();
      m_C = SEQUENCE.GetData()( middle_index)->GetAtom( GetAtomTypes().C).GetCoordinates();

      // iterate through the AAs
      for
      (
        AASequence::const_iterator aa_itr( SEQUENCE.GetData().Begin()),
          aa_itr_end( SEQUENCE.GetData().End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // get the phi/psi angles
        m_Angles.PushBack( SEQUENCE.CalculatePhiPsi( ( *aa_itr)->GetSeqID()));
      }
    }

  } // namespace biol

} // namespace bcl
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
#include "biol/bcl_biol_aa_side_chain_factory.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_classes.h"
#include "biol/bcl_biol_aa_complete.h"
#include "io/bcl_io_file.h"
#include "quality/bcl_quality_rmsd.h"
#include "score/bcl_score.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AASideChainFactory::s_Instance
    (
      GetObjectInstances().AddInstance( new AASideChainFactory())
    );

    //! @brief filename of default atom coordinate table
    const std::string &AASideChainFactory::GetDefaultTableFileName()
    {
      // initialize static table filename
      static const std::string s_filename( "idealized_side_chain_coords.table");

      // end
      return s_filename;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //!@brief default constructor
    AASideChainFactory::AASideChainFactory() :
      m_IncludeHydrogen(),
      m_IncludeBackBone(),
      m_TableFilename()
    {
    }

    //! @brief default constructor
    //! @param INCLUDE_HYDROGEN = true, if hydrogen atoms will be included in attached side chains
    //! @param INCLUDE_BACKBONE = use the backbone to align sidechain for first chi angle (important for proline)
    //! @param TABLE_FILE_NAME filename to table which contains ideal side chain conformations
    AASideChainFactory::AASideChainFactory
    (
      bool INCLUDE_HYDROGEN,
      const bool INCLUDE_BACKBONE,
      const std::string &TABLE_FILE_NAME //= GetDefaultTableFileName()
     ) :
      m_IncludeHydrogen( INCLUDE_HYDROGEN),
      m_IncludeBackBone( INCLUDE_BACKBONE),
      m_TableFilename( TABLE_FILE_NAME)
    {
      InitializeTable();
    }

    //! @brief virtual copy constructor
    AASideChainFactory *AASideChainFactory::Clone() const
    {
      return new AASideChainFactory( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AASideChainFactory::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add side chains to AABase
    //! @param AA_BASE object to which side chain atoms should be added
    void AASideChainFactory::AddSideChain( AABase &AA_BASE) const
    {
      // Generate idealized AA sidechain from AABase's type
      util::ShPtrVector< Atom> side_chain_atoms( GenerateAtoms( AA_BASE.GetType()));

      // determine the transformation matrix to put sidechains in correct place
      math::TransformationMatrix3D transformation;

      if( m_IncludeBackBone)
      {
        transformation = GetQuadTransformationMatrix( AA_BASE, side_chain_atoms);
      }
      else
      {
        transformation = GetTransformationMatrix( AA_BASE, side_chain_atoms);
      }

      // check if transformation is defined
      if( !transformation.IsDefined())
      {
        BCL_MessageStd
        (
          "undefined transformation matrix - cannot attach sidechain atoms with defined coordinates"
        );

        // iterate over all sidechain atoms and render coordiantes undefined (nan)
        for
        (
          util::ShPtrVector< Atom>::iterator
            this_itr( side_chain_atoms.Begin()), end_itr( side_chain_atoms.End());
          this_itr != end_itr; ++this_itr
        )
        {
          ( *this_itr)->SetCoordinates( linal::Vector3D( util::GetUndefined< double>()));
        }

        // set the atoms
        AA_BASE.SetAtoms( side_chain_atoms);

        // end
        return;
      }

      // transform the atoms
      TransformAtoms( transformation, side_chain_atoms);

      if( m_IncludeBackBone)
      {
        // storage for new side chain atoms
        util::ShPtrVector< Atom> side_chain_atoms_new;

        // iterate over the atoms in side_chain_atoms
        for
        (
          util::ShPtrVector< Atom>::iterator
            this_itr( side_chain_atoms.Begin()), end_itr( side_chain_atoms.End());
          this_itr != end_itr; ++this_itr
        )
        {
          // if not a backbone atom type
          if( !( *this_itr)->GetType()->IsBackBone())
          {
            // add it to side chain atoms
            side_chain_atoms_new.PushBack( *this_itr);
          }
        }
        // update the side chain atoms
        side_chain_atoms = side_chain_atoms_new;
      }

      // set the atoms
      AA_BASE.SetAtoms( side_chain_atoms);
    }

    // takes a Protein model and returns a new protein model with side chains attached to new amino acids, sses and chains
    util::ShPtr< assemble::ProteinModel> AASideChainFactory::ProteinModelWithSideChains
    (
      const assemble::ProteinModel &THIS_MODEL
    ) const
    {
      // make ShPtrVector for the new chains
      util::ShPtrVector< assemble::Chain> new_chains;

      // iterate over chains in THIS_MODEL, add side chains to the AAs in it, and push it back into new_chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( THIS_MODEL.GetChains().Begin()), chain_end_itr( THIS_MODEL.GetChains().End());
        chain_itr != chain_end_itr;
        ++chain_itr
      )
      {
        new_chains.PushBack( NewChain( **chain_itr));
      }

      // construct a new ProteinModel and return it
      util::ShPtr< assemble::ProteinModel> new_model( new assemble::ProteinModel( new_chains));
      return new_model;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AASideChainFactory::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_IncludeHydrogen, ISTREAM);
      io::Serialize::Read( m_IncludeBackBone, ISTREAM);
      io::Serialize::Read( m_Table, ISTREAM);
      io::Serialize::Read( m_TableFilename, ISTREAM);

      //return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &AASideChainFactory::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_IncludeHydrogen, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IncludeBackBone, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Table, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TableFilename, OSTREAM, INDENT);

      //return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief initialize table with idealized side chain coordinates
    void AASideChainFactory::InitializeTable()
    {
      io::IFStream read;
      io::File::MustOpenIFStream( read, score::Score::AddHistogramPath( m_TableFilename));
      m_Table.ReadFormatted( read);
      io::File::CloseClearFStream( read);
    }

    //! @brief takes a chain and returns a new chain with side chains attached to new amino acids
    //! @param THIS_CHAIN the chain which contains AAs with no side chains
    //! @return ShPtr to new chain with side chains on the AAs
    util::ShPtr< assemble::Chain> AASideChainFactory::NewChain( const assemble::Chain &THIS_CHAIN) const
    {
      util::ShPtrVector< assemble::SSE> new_sses;
      // iterate over chain and Get SSEs.
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( THIS_CHAIN.GetData().Begin()), sse_end_itr( THIS_CHAIN.GetData().End());
        sse_itr != sse_end_itr; ++sse_itr
      )
      {
        new_sses.PushBack( NewSSE( **sse_itr));
      }
      util::ShPtr< assemble::Chain> new_chain( new assemble::Chain( THIS_CHAIN.GetSequence(), new_sses));
      return new_chain;
    }

    //! @brief similar to assemble::SSE::HardCopy, but adds SideChain
    //! @param THIS_SSE the SSE to be copied and have side chains added
    //! @return ShPtr to new SSE with side chains
    util::ShPtr< assemble::SSE> AASideChainFactory::NewSSE( const assemble::SSE &THIS_SSE) const
    {
      // is the SSE considered a coil?
      //const bool sse_is_coil( THIS_SSE.GetType() == biol::GetSSTypes().COIL);

      // storage for the new amino acids with side chains attached
      util::ShPtrVector< AABase> new_aabase_vector;
      for
      (
        AASequence::const_iterator
          aa_itr( THIS_SSE.Begin()), aa_end_itr( THIS_SSE.End());
        aa_itr != aa_end_itr; ++aa_itr
      )
      {
         // create a new amino acid from this amino acid with side chains
         new_aabase_vector.PushBack( NewAAWithSideChain( **aa_itr));
      }

      AASequence new_sequence( new_aabase_vector, THIS_SSE.GetChainID(), THIS_SSE.GetFastaHeader());
      util::ShPtr< assemble::SSE> new_sse( new assemble::SSE( new_sequence, THIS_SSE.GetType()));
      return new_sse;
    }

    //! @brief adds a side chain to an individual amino acid
    //! @param THIS_AA the amino acid that should have a side chain attached to it
    //! @return ShPtr to an AAComplete with sidechain attached
    util::ShPtr< AABase> AASideChainFactory::NewAAWithSideChain( const AABase &THIS_AA) const
    {
      // make sure that THIS_AA is at least AABackBone or AAComplete
      BCL_Assert
      (
        THIS_AA.GetAAClass() == GetAAClasses().e_AABackBone ||
        THIS_AA.GetAAClass() == GetAAClasses().e_AAComplete,
        "AA Could not have side chain added because it was of type: " +
        util::Format()( THIS_AA.GetAAClass())
      );

      // create new ShPtr to new AAComplete, and add a sidechain to it
      util::ShPtr< AABase> new_aa_shptr( new AAComplete( THIS_AA));
      AddSideChain( *new_aa_shptr);
      return new_aa_shptr;
    }

    //! @brief pull the atoms for a given aa type from the table with their idealized coordinates
    //! @brief basically this is the idealized amino acid copy in space
    //! @param AATYPE the AAType considered
    util::ShPtrVector< Atom> AASideChainFactory::GenerateAtoms( const AAType &AATYPE) const
    {
      // create storage for atoms
      util::ShPtrVector< Atom> atoms;

      // grab the row corresponding to the requested AAType
      const storage::Row< double> &aa_row( m_Table[ AATYPE->GetThreeLetterCode()]);

      // get set of atom types for the AATYPE
      const storage::Set< AtomType> this_atomtypes( AATYPE->GetAllowedAtomTypes());

      // iterate over all the atom types
      for
      (
        storage::Set< AtomType>::const_iterator itr( this_atomtypes.Begin()), itr_end( this_atomtypes.End());
        itr != itr_end;
        ++itr
      )
      {
        // if hydrogens are excluded and the atom type is a hydrogen and the aatype is glycine, still include the HA
        if( AATYPE == GetAATypes().GLY && !m_IncludeHydrogen && ( *itr) == GetAtomTypes().HA)
        {
          //Add the HA for glycine
          Atom this_atom( GenerateAtomWithCoords( *itr, aa_row));

          // insert atom into m_RotamAtoms
          atoms.PushBack( util::ShPtr< Atom>( this_atom.Clone()));

          // go to next atom type
          continue;
        }

        // if hydrogens are not to be included
        if( !m_IncludeHydrogen && ( *itr)->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
        {
          continue;
        }

        // exclude backbone if desired
        if( !m_IncludeBackBone && ( *itr)->IsBackBone())
        {
          continue;
        }

        //generate atom for that atomtype using default coordinates from the table
        Atom this_atom( GenerateAtomWithCoords( *itr, aa_row));

        // insert atom into m_RotamAtoms
        atoms.PushBack( util::ShPtr< Atom>( this_atom.Clone()));
      }

      return atoms;
    }

    //! @brief transforms ShPtrVector< Atom> according to transformation matrix
    //! @param TRANSFORMATION matrix used to transform atoms
    //! @param ATOMS atoms to be transformed
    void AASideChainFactory::TransformAtoms
    (
      const math::TransformationMatrix3D &TRANSFORMATION, util::ShPtrVector< Atom> &ATOMS
    )
    {
      // iterate over the atoms in ATOMS transforming each atom
      for( util::ShPtrVector< Atom>::iterator itr( ATOMS.Begin()), end_itr( ATOMS.End()); itr != end_itr; ++itr)
      {
        ( *itr)->Transform( TRANSFORMATION);
      }
    }

    //! @brief get transformation matrix to superimpose CA-CB vector of amino acid with the CA-CB vector of the atoms
    //! @param AMINO_ACID the amino acid that is to be used
    //! @param ATOMS the atoms generated from the table
    //! @return transformationmatrix the transformation matrix used to transform the atoms of the rotamer
    math::TransformationMatrix3D AASideChainFactory::GetTransformationMatrix
    (
      const AABase &AMINO_ACID,
      const util::ShPtrVector< Atom> &ATOMS
    ) const
    {
      // make SiPtrVector of resultant for rotamer Ca-Cb Vector3D
      util::SiPtrVector< linal::Vector3D> rotam_siptr_vec;
      const linal::Vector3D &aa_ca_coords( AMINO_ACID.GetAtom( GetAtomTypes().CA).GetCoordinates());
      const linal::Vector3D &aa_cb_coords( AMINO_ACID.GetAtom( GetAtomTypes().CB).GetCoordinates());
      linal::Vector3D rotam_ca_coords;
      linal::Vector3D rotam_cb_coords;

      // iterate over atoms contained in ATOMS (rotamer atoms)
      for( util::ShPtrVector< Atom>::const_iterator itr( ATOMS.Begin()), end_itr( ATOMS.End()); itr != end_itr; ++itr)
      {
        //find the CA and get the coordinates
        if( ( *itr)->GetType() == GetAtomTypes().CA)
        {
          rotam_ca_coords = ( ( *itr)->GetCoordinates());
        }

        //find the CB and get the coordinates
        if( ( *itr)->GetType() == GetAtomTypes().CB)
        {
          rotam_cb_coords = ( ( *itr)->GetCoordinates());
        }
      }
      // define the transformation matrix
      math::TransformationMatrix3D matrix;

      // apply transformation to matrix with the negative of the rotam_ca_coords so that it translates to the origin
      matrix( -rotam_ca_coords);

      // determine the projection angle for Ca-Cb vector alignment
      const double projangle( linal::ProjAngle( rotam_ca_coords, rotam_cb_coords, aa_ca_coords, aa_cb_coords));

      // instantiate the rotation matrix with the axis of rotation (via cross product) and the projection angle
      const math::RotationMatrix3D rotationmatrix
      (
        linal::CrossProduct( aa_ca_coords, aa_cb_coords, rotam_ca_coords, rotam_cb_coords), projangle
      );

      // apply transformation matrix with the rotation matrix from above
      matrix( rotationmatrix);

      // add the aa_ca_coord vector to the transformation matrix so that the Ca-Cb is properly aligned
      // note that this math:Vector3D is implicitly converted to math::TransformationMatrix3D
      matrix( aa_ca_coords);

      // return the matrix that will apply the transformation to a set of un-transformed atoms
      return matrix;
    }

    //! @brief get transformation matrix to superimpose C-CA-CB-N vector of amino acid with the C-CA-CB-N vector of the atoms
    //! @param AMINO_ACID the amino acid that is to be used for superimposition
    //! @param ATOMS the atoms generated from the table
    //! @return transformation matrix the transformation matrix used to transform the atoms of the rotamer
    math::TransformationMatrix3D AASideChainFactory::GetQuadTransformationMatrix
    (
      const AABase &AMINO_ACID,
      const util::ShPtrVector< Atom> &ATOMS
    ) const
    {
      // initialize static set of desired atom types
      static const storage::Set< AtomType> desired_atoms
      (
        storage::Set< AtomType>::Create
        (
          GetAtomTypes().CA, GetAtomTypes().CB, GetAtomTypes().N, GetAtomTypes().C
        )
      );

      // make SiPtrVector of backbone atoms in given amino acid
      util::SiPtrVector< const linal::Vector3D> aa_coords( AMINO_ACID.GetAtomCoordinates( desired_atoms));

      // make SiPtrVector of backbone atoms in rotamer atoms
      util::SiPtrVector< const linal::Vector3D> rotam_siptr_vec;

      // iterate over backbone atom types
      for
      (
        storage::Set< AtomType>::const_iterator
          aatype_itr( desired_atoms.Begin()), aatype_itr_end( desired_atoms.End());
        aatype_itr != aatype_itr_end; ++aatype_itr
      )
      {
        // iterate over atoms contained in ATOMS (rotamer atoms)
        for( util::ShPtrVector< Atom>::const_iterator itr( ATOMS.Begin()), itr_end( ATOMS.End()); itr != itr_end; ++itr)
        {
          if( *aatype_itr == ( *itr)->GetType())
          {
            rotam_siptr_vec.PushBack( util::ToSiPtr( ( *itr)->GetCoordinates()));
            break;
          }
        }
      }
      // superimpose coordiantes and return the superimposition matrix
      return quality::RMSD::SuperimposeCoordinates( aa_coords, rotam_siptr_vec);
    }

    //! @brief a helper function which generates atoms from given atom type with its ideal coordinates from table
    //! @param ATOM_TYPE atom type to be search table for
    //! @param ROW searches the row (predetermined AA type) for the atom name of choice
    //! @return biol::Atom with coordinates from table
    Atom AASideChainFactory::GenerateAtomWithCoords
    (
      const AtomType &ATOM_TYPE,
      const storage::Row< double> &ROW
    )
    {
      // initialize atom names
      const std::string atom_x_name( ATOM_TYPE.GetName() + "_x");
      const std::string atom_y_name( ATOM_TYPE.GetName() + "_y");
      const std::string atom_z_name( ATOM_TYPE.GetName() + "_z");

      //get the coordinates based on atom name given
      const linal::Vector3D coord
      (
        ROW[ atom_x_name],
        ROW[ atom_y_name],
        ROW[ atom_z_name]
      );

      // coordinates should be defined
      BCL_Assert( coord.IsDefined(), "Some coordinates for this atom type are not defined!");

      // construct and return atom
      return Atom( coord, ATOM_TYPE);
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_aa_type_data.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "biol/bcl_biol_atom_group_types.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_molecular_configuration_shared.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "linal/bcl_linal_matrix.h"
#include "math/bcl_math_sum_function.h"
#include "sdf/bcl_sdf_bond_info.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief PropertyType as string
    //! @param PROPERTY_TYPE the PropertyType
    //! @return the string for the PropertyType
    const std::string &AATypeData::GetPropertyDescriptor( const PropertyType &PROPERTY_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "AA_NaturalPrevalence",
        "AA_StericalParameter",
        "AA_Polarizability",
        "AA_Volume",
        "AA_Hydrophobicity",
        "AA_IsoelectricPoint",
        "AA_Mass",
        "AA_Charge",
        "AA_pK_EMBOSS",
        "AA_pK_DTASelect",
        "AA_pk_Solomon",
        "AA_pK_Sillero",
        "AA_pK_Rodwell",
        "AA_pK_Patrickios",
        "AA_pK_Wikipedia",
        "AA_pK_Lehninger",
        "AA_pK_Grimsely",
        "AA_pK_Bjellqvist",
        "AA_pK_ProMoST",
        "AA_pK_Bjellqvist_NTerm",
        "AA_pK_Bjellqvist_CTerm",
        "AA_pK_ProMoST_NTerm",
        "AA_pK_ProMoST_CTerm",
        "AA_pK_Carey_NTerm",
        "AA_pK_Carey_CTerm",
        "AA_HelixProbability",
        "AA_StrandProbability",
        "AA_FreeEnergyHelix",
        "AA_FreeEnergyStrand",
        "AA_FreeEnergyCoil",
        "AA_TransferFreeEnergyWhimleyWhite",
        "AA_TransferFreeEnergyEngelmanSeitzGoldman",
        "AA_TransferFreeEnergyKyteDoolittle",
        "AA_TransferFreeEnergyEisenberg",
        "AA_TransferFreeEnergyHoppWoods",
        "AA_TransferFreeEnergyGuy",
        "AA_TransferFreeEnergyJanin",
        "AA_TransferFreeEnergyPuntaMaritan1D",
        "AA_TransferFreeEnergyPuntaMaritan3D",
        "AA_FreeEnergyCore",
        "AA_FreeEnergyTransition",
        "AA_FreeEnergySolution",
        "AA_FreeEnergyCoreHelix",
        "AA_FreeEnergyTransitionHelix",
        "AA_FreeEnergySolutionHelix",
        "AA_FreeEnergyCoreStrand",
        "AA_FreeEnergyTransitionStrand",
        "AA_FreeEnergySolutionStrand",
        "AA_FreeEnergyCoreCoil",
        "AA_FreeEnergyTransitionCoil",
        "AA_FreeEnergySolutionCoil",
        "AA_FreeEnergyCorePore",
        "AA_FreeEnergyCoreMembrane",
        "AA_MembraneStrandOrientationHydrophobicity",
        "AA_FreeEnergyExtracellularBlastBB",
        "AA_FreeEnergyExtracellularTypeBB",
        "AA_SASA",
        "AA_SideChainGirth",
        "AA_SideChainPolarizability",
        "AA_TopologicalPolarSurfaceArea",
        "AA_VdwSurfaceArea",
        "AA_HAcceptors",
        "AA_HDonors",
        "AA_Aromatic",
        GetStaticClassName< PropertyType>()
      };

      return s_descriptors[ PROPERTY_TYPE];
    }

    //! @brief helper function to create a molecular configuration given the component atom types
    //! @param ATOM_TYPES atom types from the AA
    //! @param ONE_LETTER_CODE one letter code; needed to process the exceptions to the general AA connectivity rules
    //! @param THREE_LETTER_CODE three letter code, used for naming the molecule
    //! @param IS_NATURAL all natural AAs must be saturated; so this serves as a check
    //! @param C_TERMINAL whether to include the C-terminal OH group in the molecule
    //! @param N_TERMINAL whether to include the N-terminal H in the molecule
    //! @param HYDROGENATE_TERMINAL_CARBOXYL whether the terminal O is hydrogenated (ignored if !C_TERMINAL)
    //! @param HYDROGENATE_TERMINAL_AMIDE whether the terminal N is hydrogenated (ignored if !N_TERMINAL)
    chemistry::FragmentConfigurationShared CreateFragment
    (
      const storage::Set< AtomType> &ATOM_TYPES,
      const char &ONE_LETTER_CODE,
      const std::string &THREE_LETTER_CODE,
      const bool &IS_NATURAL,
      const bool &C_TERMINAL,
      const bool &N_TERMINAL,
      const bool &HYDROGENATE_TERMINAL_CARBOXYL,
      const bool &HYDROGENATE_TERMINAL_AMIDE
    )
    {
      // create a vector with all component atom types
      storage::Vector< AtomType> atom_type_vec( ATOM_TYPES.Begin(), ATOM_TYPES.End());
      storage::Vector< int> charge_vec( atom_type_vec.GetSize(), 0);
      if( N_TERMINAL && !ATOM_TYPES.IsEmpty())
      {
        atom_type_vec.PushBack( GetAtomTypes().H1);
        charge_vec.PushBack( 0);
        const size_t h_pos( atom_type_vec.Find( GetAtomTypes().H));
        if( h_pos < atom_type_vec.GetSize())
        {
          atom_type_vec.RemoveElements( h_pos, 1);
          charge_vec.RemoveElements( h_pos, 1);
          atom_type_vec.PushBack( GetAtomTypes().H2);
          charge_vec.PushBack( 0);
          if( HYDROGENATE_TERMINAL_AMIDE)
          {
            atom_type_vec.PushBack( GetAtomTypes().H3);
            charge_vec.PushBack( 0);
            charge_vec( atom_type_vec.Find( GetAtomTypes().N)) = 1;
          }
        }
        else if( HYDROGENATE_TERMINAL_AMIDE)
        {
          // terminal proline-like residue
          atom_type_vec.PushBack( GetAtomTypes().H2);
          charge_vec.PushBack( 0);
          charge_vec( atom_type_vec.Find( GetAtomTypes().N)) = 1;
        }
      }
      if( C_TERMINAL && !ATOM_TYPES.IsEmpty())
      {
        atom_type_vec.PushBack( GetAtomTypes().OXT);
        charge_vec.PushBack( 0);
        if( HYDROGENATE_TERMINAL_CARBOXYL)
        {
          atom_type_vec.PushBack( GetAtomTypes().HXT);
          charge_vec.PushBack( 0);
        }
        else
        {
          charge_vec.LastElement() = -1;
        }
      }

      // create a vector with just the IDs from each of those atom types
      const size_t n_atom_types( atom_type_vec.GetSize());

      // create a matrix to hold the bond types
      linal::Matrix< size_t> bond_type_matrix( n_atom_types, n_atom_types, size_t( 0));

      // track number of electrons in bonds for each atom
      storage::Vector< size_t> bond_counts( n_atom_types, size_t( 0));
      storage::Vector< size_t> is_side_chain( n_atom_types, size_t( 0));
      // handle mandatory connections
      for( size_t atom_type_index_a( 0); atom_type_index_a < n_atom_types; ++atom_type_index_a)
      {
        if( atom_type_vec( atom_type_index_a)->IsSideChain())
        {
          is_side_chain( atom_type_index_a) = 1;
        }
        for( size_t atom_type_index_b( atom_type_index_a + 1); atom_type_index_b < n_atom_types; ++atom_type_index_b)
        {
          if( atom_type_vec( atom_type_index_a)->GetConnections().Contains( atom_type_vec( atom_type_index_b)))
          {
            bond_type_matrix( atom_type_index_a, atom_type_index_b) = 1;
            bond_type_matrix( atom_type_index_b, atom_type_index_a) = 1;
            ++bond_counts( atom_type_index_a);
            ++bond_counts( atom_type_index_b);
          }
        }
      }

      // handle the two AAs with non-standard bonds
      if( ONE_LETTER_CODE == 'P')
      {
        // proline, add connection between N and CD
        const size_t atom_type_index_a( atom_type_vec.Find( GetAtomTypes().N));
        const size_t atom_type_index_b( atom_type_vec.Find( GetAtomTypes().CD));
        bond_type_matrix( atom_type_index_a, atom_type_index_b) = 1;
        bond_type_matrix( atom_type_index_b, atom_type_index_a) = 1;
        ++bond_counts( atom_type_index_a);
        ++bond_counts( atom_type_index_b);
      }
      else if( ONE_LETTER_CODE == 'H')
      {
        // histidine, need to set a double bond between CD2 and CG
        const size_t atom_type_index_a( atom_type_vec.Find( GetAtomTypes().CG));
        const size_t atom_type_index_b( atom_type_vec.Find( GetAtomTypes().CD2));
        bond_type_matrix( atom_type_index_a, atom_type_index_b) = 2;
        bond_type_matrix( atom_type_index_b, atom_type_index_a) = 2;
        ++bond_counts( atom_type_index_a);
        ++bond_counts( atom_type_index_b);
      }

      // add double bonds to all necessary atoms
      for( size_t atom_type_index_a( 0); atom_type_index_a < n_atom_types; ++atom_type_index_a)
      {
        // skip atoms that have 4 bonds already
        if( bond_counts( atom_type_index_a) == size_t( 4))
        {
          continue;
        }
        for( size_t atom_type_index_b( atom_type_index_a + 1); atom_type_index_b < n_atom_types; ++atom_type_index_b)
        {
          // skip unbonded atoms
          if( !bond_type_matrix( atom_type_index_a, atom_type_index_b))
          {
            continue;
          }

          // skip atoms that have 4 bonds already
          if( bond_counts( atom_type_index_b) == size_t( 4))
          {
            continue;
          }

          // insert the double bond if it would be made to this atom
          if( atom_type_vec( atom_type_index_a)->GetDoubleBondConnections().Contains( atom_type_vec( atom_type_index_b)))
          {
            bond_type_matrix( atom_type_index_a, atom_type_index_b) = 2;
            bond_type_matrix( atom_type_index_b, atom_type_index_a) = 2;
            ++bond_counts( atom_type_index_a);
            ++bond_counts( atom_type_index_b);
          }
        }
      }

      if( IS_NATURAL)
      {
        // fix charges on side chain atoms
        for( size_t atom_type_index( 0); atom_type_index < n_atom_types; ++atom_type_index)
        {
          // only handle side chain atoms
          if( !is_side_chain( atom_type_index))
          {
            continue;
          }

          // ignore elements that are never ionized in the natural amino acids
          const AtomType &atom_type( atom_type_vec( atom_type_index));
          const chemistry::ElementType &element( atom_type->GetElementType());
          size_t main_group( element->GetMainGroup());

          // handle nitrogen group
          if( main_group == chemistry::GetElementTypes().e_Nitrogen->GetMainGroup())
          {
            if( bond_counts( atom_type_index) == size_t( 4))
            {
              charge_vec( atom_type_index) = 1;
            }
          }
          else if( main_group == chemistry::GetElementTypes().e_Oxygen->GetMainGroup())
          {
            // handle oxygen group (O, S, Se)
            if( bond_counts( atom_type_index) == size_t( 1))
            {
              charge_vec( atom_type_index) = -1;
            }
          }
          // other groups are never ionized in the natural amino acids
        }
      }

      storage::Vector< sdf::BondInfo> bond_info;
      storage::Vector< sdf::AtomInfo> atom_info( n_atom_types);
      std::string biol_atom_types;
      for( size_t atom_type_index_a( 0); atom_type_index_a < n_atom_types; ++atom_type_index_a)
      {
        atom_info( atom_type_index_a)
          = sdf::AtomInfo
            (
              chemistry::GetAtomTypes().GetAtomType
              (
                atom_type_vec( atom_type_index_a)->GetElementType(),
                charge_vec( atom_type_index_a)
              ),
              chemistry::e_NonChiral
            );
        biol_atom_types += atom_type_vec( atom_type_index_a).GetName() + " ";
        for( size_t atom_type_index_b( atom_type_index_a + 1); atom_type_index_b < n_atom_types; ++atom_type_index_b)
        {
          const size_t bond_order( bond_type_matrix( atom_type_index_a, atom_type_index_b));
          if( bond_order)
          {
            bond_info.PushBack
            (
              sdf::BondInfo
              (
                atom_type_index_a,
                atom_type_index_b,
                bond_order == 2
                ? chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond
                : chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond
              )
            );
          }
        }
      }
      for( size_t atom_type_index_a( 0); atom_type_index_a < n_atom_types; ++atom_type_index_a)
      {
        if( bond_counts( atom_type_index_a) == size_t( 0))
        {
          BCL_MessageCrt
          (
            "Atom type: " + atom_type_vec( atom_type_index_a).GetName()
            + " has no connected atoms in " + THREE_LETTER_CODE
          );
        }
        if( atom_type_vec( atom_type_index_a)->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
        {
          BCL_Assert
          (
            bond_counts( atom_type_index_a) == size_t( 1),
            "H should only have 1 bond, but has " + util::Format()( bond_counts( atom_type_index_a))
            + " in " + THREE_LETTER_CODE
          );
        }
        else
        {
          BCL_Assert
          (
            bond_counts( atom_type_index_a) < size_t( 5),
            "Atom type: " + atom_type_vec( atom_type_index_a).GetName()
            + " has " + util::Format()( bond_counts( atom_type_index_a)) + " connected atoms in " + THREE_LETTER_CODE
          );
        }
      }

      // set chirality of CA and CB

      // set CA chirality
      if( ATOM_TYPES.Contains( GetAtomTypes().CA))
      {
        chemistry::ChiralityEnum ca_chirality( chemistry::e_SChirality);
        if( ONE_LETTER_CODE == 'C' || THREE_LETTER_CODE == "DAL" || THREE_LETTER_CODE == "DPN")
        {
          // Cysteine is the sole natural AA with R chirality
          ca_chirality = chemistry::e_RChirality;
        }
        atom_info( atom_type_vec.Find( GetAtomTypes().CA)).SetChirality( ca_chirality);
      }

      // handle CB chirality
      if( ATOM_TYPES.Contains( GetAtomTypes().CB))
      {
        chemistry::ChiralityEnum cb_chirality( chemistry::e_NonChiral);
        if( ONE_LETTER_CODE == 'T')
        {
          // threonine has R-chirality on its CB atom
          cb_chirality = chemistry::e_RChirality;
        }
        else if( ONE_LETTER_CODE == 'I')
        {
          // isoleucine has S-chirality on its CB atom
          cb_chirality = chemistry::e_SChirality;
        }
        atom_info( atom_type_vec.Find( GetAtomTypes().CB)).SetChirality( cb_chirality);
      }

      // proline has variable N-chirality
      if( ONE_LETTER_CODE == 'P')
      {
        atom_info( atom_type_vec.Find( GetAtomTypes().N)).SetChirality( chemistry::e_UnknownChirality);
      }

      chemistry::AtomVector< chemistry::AtomComplete> atoms( atom_info, bond_info);
      // determine atom and bond types
      chemistry::AtomsCompleteStandardizer( atoms, THREE_LETTER_CODE, false);
      const size_t atoms_old( atoms.GetSize());
      chemistry::HydrogensHandler::SaturatePartial( atoms, is_side_chain);
      const size_t atoms_new( atoms.GetSize());
      BCL_Assert
      (
        atoms_new == atoms_old || !IS_NATURAL,
        util::Format()( atoms_new - atoms_old) + " H should be added to " + THREE_LETTER_CODE
      );

      // create a vector with indices of all side chain atoms, and a separate vector with all non-side-chain-atom ids
      storage::Vector< size_t> side_chain_atom_ids, nonside_chain_atom_ids;
      for( size_t atom_type_index( 0); atom_type_index < n_atom_types; ++atom_type_index)
      {
        if( atom_type_vec( atom_type_index)->IsSideChain())
        {
          side_chain_atom_ids.PushBack( atom_type_index);
        }
        else
        {
          nonside_chain_atom_ids.PushBack( atom_type_index);
        }
      }
      // add HUnk atom type for all missing H
      for( size_t added_h( atoms_old); added_h < atoms_new; ++added_h)
      {
        biol_atom_types += "HUnk ";
        side_chain_atom_ids.PushBack( added_h);
      }
      chemistry::FragmentComplete fragment_complete
      (
        atoms,
        THREE_LETTER_CODE,
        storage::Map< std::string, std::string>::Create
        (
          std::pair< std::string, std::string>( "BiolAtomTypes", biol_atom_types)
        )
      );

      // calculate all conformation & neighbor AA independent descriptors and store them on the molecule
      static storage::Vector< descriptor::CheminfoProperty> s_descriptors
      (
        storage::Vector< descriptor::CheminfoProperty>::Create
        (
          descriptor::GetCheminfoProperties().calc_NAtoms,
          descriptor::GetCheminfoProperties().calc_NStereo,
          descriptor::GetCheminfoProperties().calc_NRotBond,
          descriptor::GetCheminfoProperties().calc_NRings,
          descriptor::GetCheminfoProperties().calc_NAromaticRings,
          descriptor::GetCheminfoProperties().calc_NConjugatedRings,
          descriptor::GetCheminfoProperties().calc_NNonConjugatedRings,
          descriptor::GetCheminfoProperties().calc_HbondAcceptor,
          descriptor::GetCheminfoProperties().calc_HbondDonor,
          descriptor::GetCheminfoProperties().calc_Polarizability,
          descriptor::GetCheminfoProperties().calc_MolWeight,
          descriptor::GetCheminfoProperties().calc_TopologicalPolarSurfaceArea,
          descriptor::GetCheminfoProperties().calc_EstCovalentSurfaceArea,
          descriptor::GetCheminfoProperties().calc_EstVdwSurfaceArea,
          descriptor::GetCheminfoProperties().calc_EstVdwSurfaceAreaCSD,
          descriptor::GetCheminfoProperties().calc_MolTotalFormalCharge
        )
      );
      // create the side chain fragment
      chemistry::AtomVector< chemistry::AtomComplete> sidechain_atoms
      (
        fragment_complete.GetAtomInfo(),
        fragment_complete.GetBondInfo()
      );
      sidechain_atoms.Reorder( side_chain_atom_ids);
      chemistry::FragmentComplete sidechain( sidechain_atoms, "");

      // create the non-side chain atoms
      chemistry::AtomVector< chemistry::AtomComplete> nonsidechain_atoms
      (
        fragment_complete.GetAtomInfo(),
        fragment_complete.GetBondInfo()
      );
      nonsidechain_atoms.Reorder( nonside_chain_atom_ids);
      chemistry::FragmentComplete backbone( nonsidechain_atoms, "");
      for
      (
        storage::Vector< descriptor::CheminfoProperty>::iterator
          itr_props( s_descriptors.Begin()), itr_props_end( s_descriptors.End());
        itr_props != itr_props_end;
        ++itr_props
      )
      {
        // store the property for each fragment complete
        fragment_complete.StoreProperty( itr_props->GetString(), ( *itr_props)->SumOverObject( fragment_complete));
        fragment_complete.StoreProperty( "SideChain" + itr_props->GetString(), ( *itr_props)->SumOverObject( sidechain));
        fragment_complete.StoreProperty( "BackBone" + itr_props->GetString(), ( *itr_props)->SumOverObject( backbone));
      }

      chemistry::FragmentConfigurationShared fragment( fragment_complete);
      if( util::GetMessenger().GetCurrentMessageLevel() == util::Message::e_Debug)
      {
        // write the AA fragment off to the screen
        fragment_complete.WriteMDL( util::GetLogger());
      }
      return fragment;
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AATypeData::s_Instance( GetObjectInstances().AddInstance( new AATypeData()));

    //! @brief construct undefined amino acid
    AATypeData::AATypeData() :
      m_ThreeLetterCode( "UND"),
      m_OneLetterCode( 'U'),
      m_IsNaturalAminoAcid( false),
      m_ParentType( "UND"),
      m_AtomTypes(),
      m_FirstSidechainAtomType()
    {
      // fill properties with undefined
      std::fill_n( m_Properties, size_t( s_NumberPropertyTypes), util::GetUndefined< double>());
    }

    //! @brief construct AATypeData with undefined properties
    //! @param THREE_LETTER_CODE Three letter code for this amino acid type
    //! @param ONE_LETTER_CODE One letter code for this amino acid type
    //! @param IS_NATURAL_AA boolean to indicate whether this amino acid type is one of 20 natural ones
    //! @param PARENT_TYPE name of parent AAType, this is the same as the AAType for natural amino acids
    AATypeData::AATypeData
    (
      const std::string &THREE_LETTER_CODE,
      const char ONE_LETTER_CODE,
      const bool IS_NATURAL_AA,
      const std::string &PARENT_TYPE
    ) :
      m_ThreeLetterCode( THREE_LETTER_CODE),
      m_OneLetterCode( ONE_LETTER_CODE),
      m_IsNaturalAminoAcid( IS_NATURAL_AA),
      m_ParentType( PARENT_TYPE),
      m_AtomTypes()
    {
      // fill properties with undefined
      std::fill_n( m_Properties, size_t( s_NumberPropertyTypes), util::GetUndefined< double>());
    }

    //! @brief construct amino acid from all its data
    AATypeData::AATypeData
    (
      const std::string &THREE_LETTER_CODE,
      const char ONE_LETTER_CODE,
      const bool IS_NATURAL_AA,
      const std::string &PARENT_TYPE,
      const storage::Set< AtomType> &ATOM_TYPES,
      const AtomType &FIRST_SIDECHAIN_ATOM_TYPE,
      const double NATURAL_PREVALENCE,
      const double STERICAL_PARAMETER,
      const double POLARIZABILITY,
      const double VOLUME,
      const double HYDROPHOBICITY,
      const double ISOELECTRIC_POINT,
      const double CHARGE,
      const double PK_EMBOSS,
      const double PK_DTASELECT,
      const double PK_SOLOMON,
      const double PK_SILLERO,
      const double PK_RODWELL,
      const double PK_PATRICKIOS,
      const double PK_WIKIPEDIA,
      const double PK_LEHNINGER,
      const double PK_GRIMSELY,
      const double PK_BJELLQVIST,
      const double PK_PROMOST,
      const double PK_BJELLQVIST_NTERM,
      const double PK_BJELLQVIST_CTERM,
      const double PK_PROMOST_NTERM,
      const double PK_PROMOST_CTERM,
      const double PK_CAREY_NTERM,
      const double PK_CAREY_CTERM,
      const double HELIX_PROBABILITY,
      const double STRAND_PROBABILITY,
      const double FREE_ENERGY_HELIX,
      const double FREE_ENERGY_STRAND,
      const double FREE_ENERGY_COIL,
      const double TRANSFER_FREE_ENERGY_WHIMLEY_WHITE,
      const double TRANSFER_FREE_ENERGY_ENGELMAN_SEITZ_GOLDMAN,
      const double TRANSFER_FREE_ENERGY_KYTE_DOOLITTLE,
      const double TRANSFER_FREE_ENERGY_EISENBERG,
      const double TRANSFER_FREE_ENERGY_HOPP_WOODS,
      const double TRANSFER_FREE_ENERGY_GUY,
      const double TRANSFER_FREE_ENERGY_JANIN,
      const double TRANSFER_FREE_ENERGY_PUNTA_MARITAN_1D,
      const double TRANSFER_FREE_ENERGY_PUNTA_MARITAN_3D,
      const double FREE_ENERGY_CORE,
      const double FREE_ENERGY_TRANSITION,
      const double FREE_ENERGY_SOLUTION,
      const double FREE_ENERGY_CORE_HELIX,
      const double FREE_ENERGY_TRANSITION_HELIX,
      const double FREE_ENERGY_SOLUTION_HELIX,
      const double FREE_ENERGY_CORE_STRAND,
      const double FREE_ENERGY_TRANSITION_STRAND,
      const double FREE_ENERGY_SOLUTION_STRAND,
      const double FREE_ENERGY_CORE_COIL,
      const double FREE_ENERGY_TRANSITION_COIL,
      const double FREE_ENERGY_SOLUTION_COIL,
      const double FREE_ENERGY_CORE_PORE,
      const double FREE_ENERGY_CORE_MEMBRANE,
      const double BETA_BARREL_HYDROPHOBICITY,
      const double FREE_ENERGY_EXT_BLAST_BB,
      const double FREE_ENERGY_EXT_TYPE_BB,
      const double SASA,
      const double SIDE_CHAIN_GIRTH,
      const double SIDE_CHAIN_POLARIZABILITY,
      const double SIDE_CHAIN_TOPOLOGICAL_POLAR_SURFACE_AREA,
      const double SIDE_CHAIN_VAN_DER_WAALS_SURFACE_AREA,
      const int    HBOND_ACCEPTORS,
      const int    HBOND_DONORS,
      const bool   IS_AROMATIC,
      const double VDW_RADIUS_CB
    ) :
      m_ThreeLetterCode( THREE_LETTER_CODE),
      m_OneLetterCode( ONE_LETTER_CODE),
      m_IsNaturalAminoAcid( IS_NATURAL_AA),
      m_ParentType( PARENT_TYPE),
      m_AtomTypes( ATOM_TYPES),
      m_FirstSidechainAtomType( FIRST_SIDECHAIN_ATOM_TYPE),
      m_DihedralAtoms()
    {
      m_Properties[ e_NaturalPrevalence]                      = NATURAL_PREVALENCE;
      m_Properties[ e_StericalParameter]                      = STERICAL_PARAMETER;
      m_Properties[ e_Polarizability]                         = POLARIZABILITY;
      m_Properties[ e_Volume]                                 = VOLUME;
      m_Properties[ e_Hydrophobicity]                         = HYDROPHOBICITY;
      m_Properties[ e_IsoelectricPoint]                       = ISOELECTRIC_POINT;
      m_Properties[ e_Mass]                                   = CalculateMass();
      m_Properties[ e_Charge]                                 = CHARGE;
      m_Properties[ e_pK_EMBOSS]                              = PK_EMBOSS;
      m_Properties[ e_pK_DTASelect]                           = PK_DTASELECT;
      m_Properties[ e_pk_Solomon]                             = PK_SOLOMON;
      m_Properties[ e_pK_Sillero]                             = PK_SILLERO;
      m_Properties[ e_pK_Rodwell]                             = PK_RODWELL;
      m_Properties[ e_pK_Patrickios]                          = PK_PATRICKIOS;
      m_Properties[ e_pK_Wikipedia]                           = PK_WIKIPEDIA;
      m_Properties[ e_pK_Lehninger]                           = PK_LEHNINGER;
      m_Properties[ e_pK_Grimsely]                            = PK_GRIMSELY;
      m_Properties[ e_pK_Bjellqvist]                          = PK_BJELLQVIST;
      m_Properties[ e_pK_ProMoST]                             = PK_PROMOST;
      m_Properties[ e_pK_Bjellqvist_NTerm]                    = PK_BJELLQVIST_NTERM;
      m_Properties[ e_pK_Bjellqvist_CTerm]                    = PK_BJELLQVIST_CTERM;
      m_Properties[ e_pK_ProMoST_NTerm]                       = PK_PROMOST_NTERM;
      m_Properties[ e_pK_ProMoST_CTerm]                       = PK_PROMOST_CTERM;
      m_Properties[ e_pK_Carey_NTerm]                         = PK_CAREY_NTERM;
      m_Properties[ e_pK_Carey_CTerm]                         = PK_CAREY_CTERM;
      m_Properties[ e_HelixProbability]                       = HELIX_PROBABILITY;
      m_Properties[ e_StrandProbability]                      = STRAND_PROBABILITY;
      m_Properties[ e_FreeEnergyHelix]                        = FREE_ENERGY_HELIX;
      m_Properties[ e_FreeEnergyStrand]                       = FREE_ENERGY_STRAND;
      m_Properties[ e_FreeEnergyCoil]                         = FREE_ENERGY_COIL;
      m_Properties[ e_TransferFreeEnergyWhimleyWhite]         = TRANSFER_FREE_ENERGY_WHIMLEY_WHITE;
      m_Properties[ e_TransferFreeEnergyEngelmanSeitzGoldman] = TRANSFER_FREE_ENERGY_ENGELMAN_SEITZ_GOLDMAN;
      m_Properties[ e_TransferFreeEnergyKyteDoolittle]        = TRANSFER_FREE_ENERGY_KYTE_DOOLITTLE;
      m_Properties[ e_TransferFreeEnergyEisenberg]            = TRANSFER_FREE_ENERGY_EISENBERG;
      m_Properties[ e_TransferFreeEnergyHoppWoods]            = TRANSFER_FREE_ENERGY_HOPP_WOODS;
      m_Properties[ e_TransferFreeEnergyGuy]                  = TRANSFER_FREE_ENERGY_GUY;
      m_Properties[ e_TransferFreeEnergyJanin]                = TRANSFER_FREE_ENERGY_JANIN;
      m_Properties[ e_TransferFreeEnergyPuntaMaritan1D]       = TRANSFER_FREE_ENERGY_PUNTA_MARITAN_1D;
      m_Properties[ e_TransferFreeEnergyPuntaMaritan3D]       = TRANSFER_FREE_ENERGY_PUNTA_MARITAN_3D;
      m_Properties[ e_FreeEnergyCore]                         = FREE_ENERGY_CORE;
      m_Properties[ e_FreeEnergyTransition]                   = FREE_ENERGY_TRANSITION;
      m_Properties[ e_FreeEnergySolution]                     = FREE_ENERGY_SOLUTION;
      m_Properties[ e_FreeEnergyCoreHelix]                    = FREE_ENERGY_CORE_HELIX;
      m_Properties[ e_FreeEnergyTransitionHelix]              = FREE_ENERGY_TRANSITION_HELIX;
      m_Properties[ e_FreeEnergySolutionHelix]                = FREE_ENERGY_SOLUTION_HELIX;
      m_Properties[ e_FreeEnergyCoreStrand]                   = FREE_ENERGY_CORE_STRAND;
      m_Properties[ e_FreeEnergyTransitionStrand]             = FREE_ENERGY_TRANSITION_STRAND;
      m_Properties[ e_FreeEnergySolutionStrand]               = FREE_ENERGY_SOLUTION_STRAND;
      m_Properties[ e_FreeEnergyCoreCoil]                     = FREE_ENERGY_CORE_COIL;
      m_Properties[ e_FreeEnergyTransitionCoil]               = FREE_ENERGY_TRANSITION_COIL;
      m_Properties[ e_FreeEnergySolutionCoil]                 = FREE_ENERGY_SOLUTION_COIL;
      m_Properties[ e_FreeEnergyCoreFacingPore]               = FREE_ENERGY_CORE_PORE;
      m_Properties[ e_FreeEnergyCoreFacingMembrane]           = FREE_ENERGY_CORE_MEMBRANE;
      m_Properties[ e_MembraneStrandOrientationHydrophobicity]= BETA_BARREL_HYDROPHOBICITY;
      m_Properties[ e_FreeEnergyExtracellularBlastBB]         = FREE_ENERGY_EXT_BLAST_BB;
      m_Properties[ e_FreeEnergyExtracellularTypeBB]          = FREE_ENERGY_EXT_TYPE_BB;
      m_Properties[ e_SASA]                                   = SASA;
      m_Properties[ e_SideChainGirth]                         = SIDE_CHAIN_GIRTH;
      m_Properties[ e_SideChainPolarizability]                = SIDE_CHAIN_POLARIZABILITY;
      m_Properties[ e_SideChainTopologicalPolarSurfaceArea]   = SIDE_CHAIN_TOPOLOGICAL_POLAR_SURFACE_AREA;
      m_Properties[ e_SideChainVanDerWaalsSurfaceArea]        = SIDE_CHAIN_VAN_DER_WAALS_SURFACE_AREA;
      m_Properties[ e_SideChainHBondAcceptors]                = HBOND_ACCEPTORS;
      m_Properties[ e_SideChainHBondDonors]                   = HBOND_DONORS;
      m_Properties[ e_Aromatic]                               = IS_AROMATIC ? 1 : 0;

      m_ChemistryTypes.Resize( GetAtomTypes().GetEnumCount());
      m_VdwRadii.Resize( m_ChemistryTypes.GetSize(), util::GetUndefinedDouble());
      // get both a fully-terminated fragment and a fully integrated fragment. Together, they have the
      // full set of atom types
      const chemistry::FragmentConfigurationShared &terminal_fragment( GetFragment( true, true));
      const chemistry::FragmentConfigurationShared &std_fragment( GetFragment( false, false));
      const storage::Vector< std::string> term_biol_types
      (
        util::SplitString( terminal_fragment.GetMDLProperty( "BiolAtomTypes"))
      );
      auto itr_biol_atom_types_str( term_biol_types.Begin());
      for( auto itr_atom( terminal_fragment.GetAtomsIterator()); itr_atom.NotAtEnd(); ++itr_atom, ++itr_biol_atom_types_str)
      {
        if( !GetAtomTypes().HaveEnumWithName( *itr_biol_atom_types_str))
        {
          continue;
        }
        AtomType atom_type( *itr_biol_atom_types_str);
        m_ChemistryTypes( atom_type.GetIndex()) = itr_atom->GetAtomType();
      }
      const storage::Vector< std::string> std_biol_types
      ( 
        util::SplitString( std_fragment.GetMDLProperty( "BiolAtomTypes"))
      );
      itr_biol_atom_types_str = std_biol_types.Begin();
      for( auto itr_atom( std_fragment.GetAtomsIterator()); itr_atom.NotAtEnd(); ++itr_atom, ++itr_biol_atom_types_str)
      {
        if( !GetAtomTypes().HaveEnumWithName( *itr_biol_atom_types_str))
        {
          continue;
        }
        AtomType atom_type( *itr_biol_atom_types_str);
        m_ChemistryTypes( atom_type.GetIndex()) = itr_atom->GetAtomType();
        m_VdwRadii( atom_type.GetIndex())
          = itr_atom->GetAtomType()->GetAtomTypeProperty( chemistry::AtomTypeData::e_VdWaalsRadiusCSD);
      }
      // backbone vdw radii; computed ignoring residues on the same SSE (except for coils) and only considering other
      // backbone and Cb atoms, also ignoring backbone O h-bonded to glycine H
      m_VdwRadii( GetAtomTypes().C)   = 1.62;
      m_VdwRadii( GetAtomTypes().CA)  = 1.68;
      m_VdwRadii( GetAtomTypes().N)   = 1.55;
      m_VdwRadii( GetAtomTypes().O)   = 0.80;
      if( VDW_RADIUS_CB)
      {
        m_VdwRadii( m_FirstSidechainAtomType) = VDW_RADIUS_CB;
      }
    }

    //! @brief virtual copy constructor
    AATypeData *AATypeData::Clone() const
    {
      return new AATypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AATypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the parent amino acid type
    //! this maps uncommon aas to their natural counterpart like Seleno-methionin to methionin
    //! @return AAType of parent for this uncommon amino acid
    const AAType &AATypeData::GetParentType() const
    {
      return *m_ParentTypePtr;
    }

    //! @brief atom types for that amino acid type without hydrogens
    //! @return a set of atom types that are part of that amino acid type without hydrogens
    storage::Set< AtomType> AATypeData::GetAllowedHeavyAtomTypes() const
    {
      // initialize set
      storage::Set< AtomType> atom_types;

      // iterate over atom types
      for
      (
        storage::Set< AtomType>::const_iterator itr( m_AtomTypes.Begin()), itr_end( m_AtomTypes.End());
        itr != itr_end; ++itr
      )
      {
        // add to set if not hydrogen
        if( ( *itr)->GetElementType() != chemistry::GetElementTypes().e_Hydrogen)
        {
          atom_types.Insert( *itr);
        }
      }

      return atom_types;
    }

    //! @brief return the pdb atom name for the given AtomType
    //! @param ATOM_TYPE atom type
    //! @return the pdb atom name with proper spacing, so that it can be written to the pdb file, empty if this amino acid type does not have this atom type
    const std::string &AATypeData::GetPDBAtomName( const AtomType &ATOM_TYPE) const
    {
      // check if the atom type is part of the amino acid type
      if( !m_AtomTypes.Contains( ATOM_TYPE))
      {
        // check if it is a terminal type
        const storage::Map< AtomType, std::string>::const_iterator term_itr
        (
          GetAtomTypes().GetTerminalExtraAtomTypes().Find( ATOM_TYPE)
        );

        // could be found, return the pdb name
        if( term_itr != GetAtomTypes().GetTerminalExtraAtomTypes().End())
        {
          return term_itr->second;
        }

        static const std::string s_undefined;
        return s_undefined;
      }

      // find entry for that atom type in map of mapped types
      const storage::Map< AtomType, std::string>::const_iterator find_itr( m_PDBAtomName.Find( ATOM_TYPE));

      // special pdb atom name defined
      if( find_itr != m_PDBAtomName.End())
      {
        return find_itr->second;
      }

      // end
      return ATOM_TYPE->AtomTypeData::GetName();
    }

    //! @brief return the associated atom type from a given atom name (IUPAC or PDB)
    //! @param ATOM_NAME the atom name defined by either the pdb or IUPAC for this aa type
    //! @return the atom type for that pdb atom name, if it is defined for that amino acid type
    AtomType AATypeData::GetAtomTypeFromAtomName( const std::string &ATOM_NAME) const
    {
      // check if AtomType can be constructed from PDB_ATOM_NAME
      const AtomType type( GetAtomTypes().TypeFromPDBAtomName( ATOM_NAME));

      // valid IUPAC name
      if( type.IsDefined())
      {
        // is pdb atom name a valid atom for this type or is atom type a terminal atom type
        if
        (
             ( m_AtomTypes.Find( type) != m_AtomTypes.End())
          || GetAtomTypes().GetTerminalExtraAtomTypes().Has( type)
        )
        {
          return type;
        }

        // wrong atom type for this aa type
        return GetAtomTypes().e_Undefined;
      }

      // iterate over pdb atom map
      for
      (
        storage::Map< AtomType, std::string>::const_iterator itr( m_PDBAtomName.Begin()), itr_end( m_PDBAtomName.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->second.find( ATOM_NAME) != std::string::npos)
        {
          return itr->first;
        }
      }

      // could still be a pdb terminal atom type name
      const AtomType term_type( GetAtomTypes().GetTerminalAtomTypeFromName( ATOM_NAME));

      // will be defined, if it is a pdb terminal atom type, or undefined otherwise
      return term_type;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief map a given atom type to a pdb atom name - is only used in the constructor of the AtomTypes class
    //! @param ATOM_TYPE the atom type
    //! @param PDB_ATOM_NAME the atom name defined by the pdb for this aa type
    //! @return true if mapping was successful (false if the atom type was already added)
    bool AATypeData::AddAtomTypePDBAtomNameMapping( const AtomType &ATOM_TYPE, const std::string &PDB_ATOM_NAME)
    {
      return m_PDBAtomName.Insert( std::pair< AtomType, std::string>( ATOM_TYPE, PDB_ATOM_NAME)).second;
    }

    //! @brief check if given ATOM_TYPE is part of the amino acid
    //! @param ATOM_TYPE type that is checked for
    //! @return true is this amino acid type has the given ATOM_TYPE
    bool AATypeData::DoesContainAtomType( const AtomType &ATOM_TYPE) const
    {
      return m_AtomTypes.Find( ATOM_TYPE) != m_AtomTypes.End();
    }

    //! @brief get properties to be used in ANN as a vector
    //! @return properties to be used in ANN as a vector
    linal::Vector< double> AATypeData::GetPropertiesForANN() const
    {
      // initialize vector of properties
      linal::Vector< double> properties_for_ANN( s_NumberPropertyTypesForANN);

      // initialize ptr to first of property vector
      double *ptr = properties_for_ANN.Begin();

      //Get 7 AA properties
      *ptr++ = m_Properties[ e_StericalParameter]; // sterical parameter
      *ptr++ = m_Properties[ e_Polarizability];    // polarizability
      *ptr++ = m_Properties[ e_Volume];            // volume
      *ptr++ = m_Properties[ e_Hydrophobicity];    // hydrophobicity
      *ptr++ = m_Properties[ e_IsoelectricPoint];  // isoelectric point
      *ptr++ = m_Properties[ e_HelixProbability];  // helix probability
      *ptr++ = m_Properties[ e_StrandProbability]; // strand probability

      // return
      return properties_for_ANN;
    }

    //! @brief get a chemistry::FragmentConfigurationalShared representation of the AA
    //! @param C_TERMINAL whether to return a representation with the C-terminal carboxyl group
    //! @param N_TERMINAL whether to return a representation with the N-terminal amino group
    //! @return a chemistry::FragmentConfigurationalShared representation of the AA
    const chemistry::FragmentConfigurationShared &AATypeData::GetFragment
    (
      const bool &C_TERMINAL,
      const bool &N_TERMINAL
    ) const
    {
      // create static representations of this AA as a fragment
      // 1. The isolated AA (fully saturated, with terminal amino/carboxyl groups)
      // 2. The AA at the N-terminus of the molecule
      // 3. The AA at the C-terminus of the molecule
      // 4. The AA with both sides having valences for a continued protein chain
      static storage::Map< std::string, storage::VectorND< 4, chemistry::FragmentConfigurationShared> > s_map;
      storage::Map< std::string, storage::VectorND< 4, chemistry::FragmentConfigurationShared> >::const_iterator
        itr( s_map.Find( m_ThreeLetterCode));
      static float s_default_ph( 7.0);
      if( itr == s_map.End())
      {
        // choose which atoms to include based on the pH
        storage::Set< AtomType> atom_types_at_ph( m_AtomTypes);

        // check for a side chain H that is ionized at the default ph
        if( m_SideChainIonizableHydrogen.IsDefined() && s_default_ph > m_Properties[ e_pK_Rodwell])
        {
          BCL_MessageVrb( "Ionizing " + m_SideChainIonizableHydrogen.GetName() + " from default " + m_ThreeLetterCode);
          atom_types_at_ph.Erase( m_SideChainIonizableHydrogen);
        }

        // determine whether to hydrogenate the terminal carboxyl
        const bool hydrogenate_term_c( s_default_ph < m_Properties[ e_pK_Carey_CTerm]);
        const bool hydrogenate_term_n( s_default_ph < m_Properties[ e_pK_Carey_NTerm]);

        // need to create the fragments
        s_map[ m_ThreeLetterCode] =
          storage::VectorND< 4, chemistry::FragmentConfigurationShared>
          (
            CreateFragment( atom_types_at_ph, m_OneLetterCode, m_ThreeLetterCode, m_IsNaturalAminoAcid, true,  true,  hydrogenate_term_c, hydrogenate_term_n),
            CreateFragment( atom_types_at_ph, m_OneLetterCode, m_ThreeLetterCode, m_IsNaturalAminoAcid, true,  false, hydrogenate_term_c, hydrogenate_term_n),
            CreateFragment( atom_types_at_ph, m_OneLetterCode, m_ThreeLetterCode, m_IsNaturalAminoAcid, false, true,  hydrogenate_term_c, hydrogenate_term_n),
            CreateFragment( atom_types_at_ph, m_OneLetterCode, m_ThreeLetterCode, m_IsNaturalAminoAcid, false, false, hydrogenate_term_c, hydrogenate_term_n)
          );
        itr = s_map.Find( m_ThreeLetterCode);
      }
      return itr->second( C_TERMINAL ? ( N_TERMINAL ? 0 : 1) : ( N_TERMINAL ? 2 : 3));
    }

    //! @brief Get the chemistry atom type of a particular atom in this aa type
    //! @param ATOM the biol::AtomType of interest
    //! @return chemistry::AtomType of that Atom
    chemistry::AtomType AATypeData::GetChemistryAtomType( const AtomType &ATOM_TYPE) const
    {
      return ATOM_TYPE.IsDefined() ? m_ChemistryTypes( ATOM_TYPE.GetIndex()) : chemistry::AtomType();
    }

    //! @brief Get the effective VdW radius of an atom type for this AA to other AAs
    //! @param ATOM the biol::AtomType of interest
    //! @return van-der-waals radius to external AAs of the given atom type
    double AATypeData::GetVdwRadiusToOtherAA( const AtomType &ATOM_TYPE) const
    {
      return ATOM_TYPE.IsDefined() ? m_VdwRadii( ATOM_TYPE.GetIndex()) : 0.0;
    }

    //! @brief gets the structure factor for this AA Type
    //! @return the structure factor for this AA Type
    util::ShPtr< math::FunctionInterfaceSerializable< restraint::SasDataParameters, double> > AATypeData::GetStructureFactor() const
    {
      // initialize structure factors
      util::ShPtr< math::SumFunction< restraint::SasDataParameters, double> > structure_factors
      (
        new math::SumFunction< restraint::SasDataParameters, double>()
      );

      // iterate over the atoms of a given amino acid
      for
      (
        storage::Set< AtomType>::const_iterator atom_itr( m_AtomTypes.Begin()), atom_itr_end( m_AtomTypes.End());
        atom_itr != atom_itr_end; ++atom_itr
      )
      {
        // Do not add hydrogen to the structure_factor function
        if( ( *atom_itr)->GetElementType() != chemistry::GetElementTypes().e_Hydrogen)
        {
          // add the structure factor
          *structure_factors += *( GetAtomGroupTypes().GetType( GetAATypes().AATypeFromThreeLetterCode( m_ThreeLetterCode), *atom_itr));
        }
      }

      // end
      return structure_factors;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AATypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ThreeLetterCode       , ISTREAM);
      io::Serialize::Read( m_OneLetterCode         , ISTREAM);
      io::Serialize::Read( m_IsNaturalAminoAcid    , ISTREAM);
      io::Serialize::Read( m_ParentType            , ISTREAM);
      io::Serialize::Read( m_AtomTypes             , ISTREAM);
      io::Serialize::Read( m_PDBAtomName           , ISTREAM);
      io::Serialize::Read( m_FirstSidechainAtomType, ISTREAM);
      io::Serialize::Read( m_DihedralAtoms         , ISTREAM);

      BCL_MessageCrt( "AAProperty reading not implemented yet");

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &AATypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ThreeLetterCode       , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_OneLetterCode         , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_IsNaturalAminoAcid    , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_ParentType            , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomTypes             , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PDBAtomName           , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FirstSidechainAtomType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DihedralAtoms         , OSTREAM, INDENT);

      BCL_MessageCrt( "AAProperty writing not implemented yet");

      // return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set the dihedral angles for this AA type
    //! @param ANGLES the dihedral angles for this AAType
    void AATypeData::SetDihedralAngles( const storage::Map< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> > &ANGLES)
    {
      m_DihedralAtoms = ANGLES;
    }

    //! @brief set the ionizable H
    //! @param ATOM_TYPE the ionizable hydrogen
    void AATypeData::SetSideChainIonizableHType( const AtomType &ATOM_TYPE)
    {
      BCL_Assert
      (
        ATOM_TYPE->GetElementType() == chemistry::GetElementTypes().e_Hydrogen,
        "Cannot set pka of non-hydrogen atom type"
      );
      m_SideChainIonizableHydrogen = ATOM_TYPE;
    }

    //! @brief get the counts of what type and how many atoms in a given amino acid
    //! @return Return the map of type and count of atoms in a given amino acid
    storage::Map< chemistry::ElementType, size_t> AATypeData::GetElementTypeCount() const
    {
      // initialize map to hold element type counts
      storage::Map< chemistry::ElementType, size_t> map;

      // iterate over the atom types
      for
      (
        storage::Set< AtomType>::const_iterator itr( m_AtomTypes.Begin()), itr_end( m_AtomTypes.End());
        itr != itr_end; ++itr
      )
      {
        // increment the counter if the map already has the element type, otherwise set it to one
        // also ensure that the atom is not part of the backbone
        const chemistry::ElementType &element_type( ( *itr)->GetElementType());

        // if this element type has been seen
        if( map.Has( element_type))
        {
          // increment the counter
          map[ element_type]++;
        }
        // first time the element has been seen
        else
        {
          // set the counter to 1
          map[ element_type] = 1;
        }
      }

      // end
      return map;
    }

    //! @brief Calculate the mass of this amino acid (not when peptide bonded!)
    //! @return the mass of the amino acid as the sum of all atom masses; non-monoisotopic
    double AATypeData::CalculateMass() const
    {
      double mass( 0.0);
      for
      (
        storage::Set< AtomType>::const_iterator itr( m_AtomTypes.Begin()), itr_end( m_AtomTypes.End());
        itr != itr_end;
        ++itr
      )
      {
        mass += ( *itr)->GetElementType()->GetProperty( chemistry::ElementTypeData::e_Mass);
      }

      return mass;
    }

    //! @brief determine the parent type of amino acid
    void AATypeData::DetermineParentType( const AATypes &AATYPES)
    {
      for
      (
        AATypes::const_iterator itr( AATYPES.Begin()), itr_end( AATYPES.End());
        itr != itr_end;
        ++itr
      )
      {
        if( m_ParentType == ( *itr)->GetThreeLetterCode())
        {
          m_ParentTypePtr = util::ToSiPtr( *itr);
          break;
        }
      }
      if( !m_ParentTypePtr.IsDefined())
      {
        m_ParentTypePtr = util::ToSiPtr( AATYPES.XXX);
      }
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_aa_types.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief construct all AATypes
    AATypes::AATypes() :
//                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          Properties                                       charge pK     pk    pk    pk     pk    pk     pk     pk    pk      pk   pk     pk    pK     pk     pk     pk    pk   SSEProb.      Free Energies           Transfer Free Energies                                                  Free Energies continued
//                                                                   3L    1L  Natural Parent                   Atoms                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   type of first SCA   Prev   Steric  Pol     Vol   Hydrop. Isoe.              EMBOSS DTA   Sol   Sil    Rod   Pat    Wik    Leh   Gri     Bje  ProM   BjeN  BjeC   ProN   ProC   CarN  CarC     H      S      H       S       C       WW      ESG     KD      EISEN   HW      GUY     JANIN   PM1D    PM3D    Core    Trans   Sol   CoreH   TransH  SolH    CoreS   TransS  SolC    CoreC   TransC  SolC  CorePore CoreMem BBHydro   ExtBl  ExtTy     SASA  Girth SCPolr, SCTPSA,  VdwSA,  HAcc HDon, Aromatic Vdw Radius of C-beta
      ALA( AddEnum( "ALANINE",               /*  0 */   AATypeData( "ALA", 'A', true,  "ALA", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().H,   GetAtomTypes().HA,  GetAtomTypes().HB1, GetAtomTypes().HB2, GetAtomTypes().HB3),                                                                                                                                                                                                                                                                                                         GetAtomTypes().CB,  0.085, 1.280, 0.050, 1.000,  0.310,  6.110,/*A*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.59, 3.55,  7.58,  3.75,  9.69, 2.34, 0.420, 0.230, -0.168,  0.097,  0.116,  0.500, -1.600,  1.800,  0.620,  0.500,  0.100, -0.300, -0.170, -0.150, -0.130,  0.075,  0.081, -0.207, -0.071, -0.136, -0.109,  0.215,  0.211, -0.044,  0.177,  0.153, -2.588,  2.327,  1.454, -0.259, -0.597, 209.020, 2.15,  5.361,  40.07,  63.35,    0,   0, false, 1.61))), //ALA
      ARG( AddEnum( "ARGININE",              /*  1 */   AATypeData( "ARG", 'R', true,  "ARG", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD,  GetAtomTypes().NE,  GetAtomTypes().CZ,  GetAtomTypes().NH1,  GetAtomTypes().NH2,  GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,   GetAtomTypes().HD2,  GetAtomTypes().HD3,   GetAtomTypes().HE,  GetAtomTypes().HH11, GetAtomTypes().HH12, GetAtomTypes().HH21, GetAtomTypes().HH22), GetAtomTypes().CB,  0.051, 2.340, 0.290, 6.130, -1.010, 10.740,/*R*/    1.0, 12.5, 12.0, 12.5, 12.0, 11.50, 11.2, 12.48, 12.40, 12.0, 12.00, 12.50, 7.50, 3.55, 11.50, 11.50,  9.04, 2.17, 0.360, 0.250,  0.049, -0.002, -0.043,  1.810, 12.300, -4.500, -2.530, -3.000,  1.910,  1.400,  0.370,  0.320,  0.444, -0.128, -0.147,  0.927, -0.144, -0.185, -0.020, -0.061, -0.002,  0.156, -0.027, -0.097, -0.624, -4.596, -4.274,  0.000,  0.000, 335.732, 7.36,  2.222,   0.00,  42.20,    0,   3, false, 1.84))), //ARG
      ASN( AddEnum( "ASPARAGINE",            /*  2 */   AATypeData( "ASN", 'N', true,  "ASN", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().OD1, GetAtomTypes().ND2, GetAtomTypes().H,   GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HD21, GetAtomTypes().HD22),                                                                                                                                                                                                                    GetAtomTypes().CB,  0.043, 1.600, 0.130, 2.950, -0.600,  6.520,/*N*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  7.22,  3.64,  8.80, 2.02, 0.210, 0.220,  0.157,  0.081, -0.182,  0.850,  4.800, -3.500, -0.780, -0.200,  0.480,  0.500,  0.180,  0.220,  0.397, -0.143, -0.114,  0.759,  0.055, -0.048,  0.048, -0.220,  0.241,  0.113, -0.243, -0.175, -2.750, -2.167, -2.367,  0.122, -0.210, 259.845, 4.35, 12.577,  63.64, 150.80,    1,   1, false, 1.78))), //ASN
      ASP( AddEnum( "ASPARTIC_ACID",         /*  3 */   AATypeData( "ASP", 'D', true,  "ASP", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().OD1, GetAtomTypes().OD2, GetAtomTypes().H,   GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HD2),                                                                                                                                                                                                                                          GetAtomTypes().CB,  0.058, 1.600, 0.110, 2.780, -0.770,  2.950,/*D*/   -1.0,  3.9,  4.4,  3.9,  4.0,  3.68,  4.2,  3.90,  3.86,  3.5,  4.05,  4.07, 7.50, 4.55,  3.57,  4.57,  9.60, 1.88, 0.250, 0.200,  0.169,  0.067, -0.180,  3.640,  9.200, -3.500, -0.900, -3.000,  0.780,  0.600,  0.370,  0.410,  0.574, -0.068, -0.237,  0.963,  0.227, -0.193,  0.237, -0.295,  0.074,  0.266, -0.176, -0.253, -4.159, -4.524, -2.153,  0.181, -0.366, 257.993, 4.00,  6.164,  43.09,  75.41,    2,   0, false, 1.77))), //ASP
      CYS( AddEnum( "CYSTEINE",              /*  4 */   AATypeData( "CYS", 'C', true,  "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().SG,  GetAtomTypes().H,   GetAtomTypes().HA,  GetAtomTypes().HB2, GetAtomTypes().HB3,  GetAtomTypes().HG),                                                                                                                                                                                                                                                                                     GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*C*/   -1.0,  8.5,  8.5,  8.3,  9.0,  8.33,  0.0,  8.18,  8.33,  6.8,  9.00,  8.28, 7.50, 3.55,  8.00,  9.00,  8.18, 1.96, 0.170, 0.410, -0.030,  0.245, -0.148, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150,  0.040,  0.369, -0.250, -0.125,  0.240, -0.285,  1.381,  1.205, -0.159, -0.228,  0.281, -0.353, -1.989,  0.378,  4.301,  0.094,  0.013, 240.500, 3.71,  5.222,  38.80,  56.60,    0,   0, false, 1.74))), //CYS
      GLN( AddEnum( "GLUTAMINE",             /*  5 */   AATypeData( "GLN", 'Q', true,  "GLN", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD,  GetAtomTypes().OE1, GetAtomTypes().NE2, GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,  GetAtomTypes().HE21, GetAtomTypes().HE22),                                                                                                                                                     GetAtomTypes().CB,  0.038, 1.560, 0.180, 3.950, -0.220,  5.650,/*Q*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  6.73,  3.57,  9.13, 2.17, 0.360, 0.250,  0.028, -0.046,  0.020,  0.770,  4.100, -3.500, -0.850, -0.200,  0.950,  0.700,  0.260,  0.030,  0.381, -0.039, -0.201,  0.797, -0.012, -0.270, -0.023, -0.136, -0.082,  0.170,  0.109, -0.056, -1.367, -2.601, -3.214,  0.042, -0.295, 286.761, 5.50,  7.999,  43.09,  94.99,    1,   1, false, 1.81))), //GLN
      GLU( AddEnum( "GLUTAMIC_ACID",         /*  6 */   AATypeData( "GLU", 'E', true,  "GLU", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD,  GetAtomTypes().OE1, GetAtomTypes().OE2, GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,  GetAtomTypes().HE2),                                                                                                                                                                           GetAtomTypes().CB,  0.069, 1.560, 0.150, 3.780, -0.640,  3.090,/*E*/   -1.0,  4.1,  4.4,  4.3,  4.5,  4.25,  4.2,  4.07,  4.25,  4.2,  4.45,  4.45, 7.70, 4.75,  4.15,  4.75,  9.67, 2.19, 0.420, 0.210,  0.013,  0.080, -0.082,  3.630,  8.200, -3.500, -0.740, -3.000,  0.830,  0.700,  0.150,  0.300,  0.553,  0.071, -0.318,  0.905,  0.090, -0.418,  0.167, -0.024, -0.039,  0.262,  0.139, -0.271, -1.607, -4.323, -2.144,  0.009, -1.768, 285.025, 5.05,  7.196,  40.07,  82.94,    2,   0, false, 1.84))), //GLU
      GLY( AddEnum( "GLYCINE",               /*  7 */   AATypeData( "GLY", 'G', true,  "GLY", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().H,  GetAtomTypes().HA2, GetAtomTypes().HA3),                                                                                                                                                                                                                                                                                                                                                                     GetAtomTypes().HA2, 0.075, 0.000, 0.000, 0.000,  0.000,  6.070,/*G*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  7.50,  3.70,  9.60, 2.34, 0.130, 0.150,  0.109,  0.081, -0.153,  1.150, -1.000, -0.400,  0.480,  0.000,  0.330, -0.300,  0.010,  0.080, -0.085, -0.019,  0.123, -0.048,  0.071,  0.408, -0.210,  0.113,  0.355, -0.105, -0.172, -0.100, -1.309,  0.544, -0.841, -0.122, -0.989, 185.154, 1.08,  0.774,   0.00,  29.86,    0,   0, false, 0.91))), //GLY
      HIS( AddEnum( "HISTIDINE",             /*  8 */   AATypeData( "HIS", 'H', true,  "HIS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().ND1, GetAtomTypes().CD2, GetAtomTypes().CE1, GetAtomTypes().NE2,  GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HD1,  GetAtomTypes().HD2,  GetAtomTypes().HE1,   GetAtomTypes().HE2),                                                                                                                                GetAtomTypes().CB,  0.025, 2.990, 0.230, 4.660,  0.130,  7.690,/*H*/    1.0,  6.5,  6.5,  6.0,  6.4,  6.00,  0.0,  6.04,  6.00,  6.6,  5.98,  6.08, 7.50, 3.55,  4.89,  6.89,  9.17, 1.82, 0.270, 0.300,  0.005,  0.081, -0.075,  0.110,  3.000, -3.200, -0.400,  0.500, -0.500,  0.100, -0.020,  0.060,  0.191, -0.110, -0.040,  0.155, -0.125,  0.001,  0.413, -0.180,  0.042,  0.061, -0.067, -0.099, -1.672,  0.538,  0.976, -0.237,  2.773, 290.040, 5.51,  9.716,  28.68,  98.73,    2,   2, true , 1.76))), //HIS
      ILE( AddEnum( "ISOLEUCINE",            /*  9 */   AATypeData( "ILE", 'I', true,  "ILE", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG1, GetAtomTypes().CG2, GetAtomTypes().CD1, GetAtomTypes().H,   GetAtomTypes().HA,   GetAtomTypes().HB,   GetAtomTypes().HG12, GetAtomTypes().HG13, GetAtomTypes().HG21, GetAtomTypes().HG22, GetAtomTypes().HG23, GetAtomTypes().HD11,  GetAtomTypes().HD12, GetAtomTypes().HD13),                                                                                                          GetAtomTypes().CB,  0.056, 4.190, 0.190, 4.000,  1.800,  6.040,/*I*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  7.48,  3.72,  9.60, 2.36, 0.300, 0.450, -0.102, -0.058,  0.204, -1.120, -3.100,  4.500,  1.380,  1.800, -1.130, -0.700, -0.280, -0.290, -0.208,  0.125,  0.157, -0.359,  0.012,  0.148,  0.144,  0.164, -0.172, -0.061,  0.170,  0.248, -3.010,  3.748,  2.775, -0.313, -1.551, 273.462, 5.11,  7.727,   0.00, 101.00,    0,   0, false, 1.93))), //ILE
      LEU( AddEnum( "LEUCINE",               /* 10 */   AATypeData( "LEU", 'L', true,  "LEU", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD1, GetAtomTypes().CD2, GetAtomTypes().H,   GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG,   GetAtomTypes().HD11, GetAtomTypes().HD12, GetAtomTypes().HD13, GetAtomTypes().HD21,  GetAtomTypes().HD22, GetAtomTypes().HD23),                                                                                                          GetAtomTypes().CB,  0.088, 2.590, 0.190, 4.000,  1.700,  6.040,/*L*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  7.46,  3.73,  9.60, 2.36, 0.390, 0.310, -0.164,  0.028,  0.191, -1.250, -2.800,  3.800,  1.060,  1.800, -1.180, -0.500, -0.280, -0.360, -0.187,  0.094,  0.151, -0.323, -0.089, -0.012, -0.042,  0.295,  0.024, -0.028,  0.202,  0.226, -2.199,  4.597,  3.560, -0.013, -1.913, 278.520, 4.70,  7.727,   0.00, 106.10,    0,   0, false, 1.80))), //LEU
      LYS( AddEnum( "LYSINE",                /* 11 */   AATypeData( "LYS", 'K', true,  "LYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD,  GetAtomTypes().CE,  GetAtomTypes().NZ,  GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,  GetAtomTypes().HD2,  GetAtomTypes().HD3,   GetAtomTypes().HE2,  GetAtomTypes().HE3,   GetAtomTypes().HZ1, GetAtomTypes().HZ2,  GetAtomTypes().HZ3),                                            GetAtomTypes().CB,  0.059, 1.890, 0.220, 4.770, -0.990,  9.990,/*K*/    1.0, 10.8, 10.0, 10.5, 10.4, 11.50, 11.2, 10.54, 10.50, 10.5, 10.00,  9.80, 7.50, 3.55, 10.00, 10.30,  8.95, 2.18, 0.320, 0.270,  0.143, -0.041, -0.078,  2.800,  8.800, -3.900, -1.500, -3.000,  1.400,  1.800,  0.320,  0.240,  0.724, -0.131, -0.220,  1.143,  0.033, -0.163,  0.402, -0.366, -0.114,  0.378, -0.085, -0.162, -0.827, -4.193, -5.659,  0.105, -0.818, 303.428, 6.57,  9.562,  27.64, 126.70,    0,   1, false, 1.82))), //LYS
      MET( AddEnum( "METHIONINE",            /* 12 */   AATypeData( "MET", 'M', true,  "MET", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().SD,  GetAtomTypes().CE,  GetAtomTypes().H,   GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,  GetAtomTypes().HE1,  GetAtomTypes().HE2,  GetAtomTypes().HE3),                                                                                                                                                      GetAtomTypes().CB,  0.021, 2.350, 0.220, 4.430,  1.230,  5.710,/*M*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.00, 3.55,  6.98,  3.68,  9.21, 2.28, 0.380, 0.320, -0.172,  0.098,  0.121, -0.670, -3.400,  1.900,  0.640,  1.300, -1.590, -0.400, -0.260, -0.190, -0.140,  0.038,  0.136, -0.305, -0.160, -0.009,  0.199,  0.471, -0.003, -0.119,  0.098,  0.170,  0.076,  3.655,  0.973, -0.219,  0.613, 291.524, 5.66,  8.892,  25.30, 102.80,    0,   0, false, 1.82))), //MET
      PHE( AddEnum( "PHENYLALANINE",         /* 13 */   AATypeData( "PHE", 'F', true,  "PHE", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD1, GetAtomTypes().CD2, GetAtomTypes().CE1, GetAtomTypes().CE2,  GetAtomTypes().CZ,   GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HD1,  GetAtomTypes().HD2,   GetAtomTypes().HE1,  GetAtomTypes().HE2,   GetAtomTypes().HZ),                                                                                      GetAtomTypes().CB,  0.039, 2.940, 0.290, 5.890,  1.790,  5.670,/*F*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  6.96,  3.98,  9.13, 1.83, 0.300, 0.380, -0.081, -0.043,  0.149, -1.710, -3.700,  2.800,  1.190,  2.500, -2.120, -0.500, -0.410, -0.220, -0.169, -0.034,  0.294, -0.260, -0.089,  0.347,  0.014,  0.092, -0.038, -0.108, -0.055,  0.381, -2.100,  4.498,  3.397, -0.038, -0.894, 311.302, 6.16, 12.426,   0.00, 120.40,    0,   0, true , 1.77))), //PHE
      PRO( AddEnum( "PROLINE",               /* 14 */   AATypeData( "PRO", 'P', true,  "PRO", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD,  GetAtomTypes().HA,  GetAtomTypes().HB2, GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,  GetAtomTypes().HD2,  GetAtomTypes().HD3),                                                                                                                                                                                                                     GetAtomTypes().CB,  0.046, 2.670, 0.000, 2.720,  0.720,  6.800,/*P*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 8.36, 3.55,  8.36,  3.61, 10.60, 1.99, 0.130, 0.340,  0.192,  0.550, -0.371,  0.140,  0.200, -1.600,  0.120,  0.000,  0.730,  0.300,  0.130,  0.150,  0.358, -0.103, -0.137,  0.253,  0.047,  0.115,  0.870,  0.398,  0.374, -0.005, -0.440, -0.463, -6.914, -0.757,  1.576,  0.049,  1.796, 235.409, 4.12,  5.505,   0.00,  79.84,    0,   0, false, 1.71))), //PRO
      SER( AddEnum( "SERINE",                /* 15 */   AATypeData( "SER", 'S', true,  "SER", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().OG,  GetAtomTypes().H,   GetAtomTypes().HA,  GetAtomTypes().HB2, GetAtomTypes().HB3,  GetAtomTypes().HG),                                                                                                                                                                                                                                                                                     GetAtomTypes().CB,  0.060, 1.310, 0.060, 1.600, -0.040,  5.700,/*S*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 6.93, 3.55,  6.86,  3.61,  9.15, 2.21, 0.200, 0.280,  0.089, -0.032, -0.048,  0.460, -0.600, -0.800, -0.180, -0.300,  0.520,  0.100,  0.050,  0.160,  0.050, -0.023, -0.023,  0.184,  0.125,  0.051, -0.037, -0.141,  0.043, -0.088, -0.058, -0.008, -0.995, -1.692, -3.127,  0.316,  1.191, 223.038, 2.66,  2.859,  20.23,  52.65,    1,   1, false, 1.69))), //SER
      THR( AddEnum( "THREONINE",             /* 16 */   AATypeData( "THR", 'T', true,  "THR", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().OG1, GetAtomTypes().CG2, GetAtomTypes().H,   GetAtomTypes().HA,  GetAtomTypes().HB,   GetAtomTypes().HG1,  GetAtomTypes().HG21, GetAtomTypes().HG22, GetAtomTypes().HG23),                                                                                                                                                                                                                    GetAtomTypes().CB,  0.055, 3.030, 0.110, 2.600,  0.260,  5.600,/*T*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 6.82, 3.55,  7.02,  3.57,  9.10, 2.09, 0.210, 0.360,  0.071, -0.085,  0.025,  0.250, -1.200, -0.700, -0.050,  0.400,  0.070,  0.200,  0.020, -0.080,  0.003,  0.009, -0.012,  0.055,  0.087,  0.062,  0.012, -0.169, -0.085, -0.010,  0.055,  0.042, -1.294,  0.076, -2.225,  0.411, -0.820, 243.554, 3.93,  4.694,  20.23,  72.23,    1,   1, false, 1.83))), //THR
      TRP( AddEnum( "TRYPTOPHAN",            /* 17 */   AATypeData( "TRP", 'W', true,  "TRP", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD1, GetAtomTypes().CD2, GetAtomTypes().NE1, GetAtomTypes().CE2,  GetAtomTypes().CE3,  GetAtomTypes().CZ2,  GetAtomTypes().CZ3,  GetAtomTypes().CH2,  GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,   GetAtomTypes().HB3,  GetAtomTypes().HD1,   GetAtomTypes().HE1, GetAtomTypes().HE3,  GetAtomTypes().HZ2, GetAtomTypes().HZ3,  GetAtomTypes().HH2),   GetAtomTypes().CB,  0.013, 3.210, 0.410, 8.080,  2.250,  5.940,/*W*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  7.11,  3.78,  9.39, 2.83, 0.320, 0.420, -0.028, -0.069,  0.113, -2.090, -1.900, -0.900,  0.810,  3.400, -0.510, -0.300, -0.150, -0.280,  0.025, -0.214,  0.299,  0.024, -0.232,  0.397,  0.047, -0.316,  0.115,  0.110, -0.108,  0.346, -1.870,  4.387,  2.867, -0.007,  0.444, 350.681, 7.38, 17.695,  15.79, 148.64,    0,   1, true , 1.80))), //TRP
      TYR( AddEnum( "TYROSINE",              /* 18 */   AATypeData( "TYR", 'Y', true,  "TYR", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD1, GetAtomTypes().CD2, GetAtomTypes().CE1, GetAtomTypes().CE2,  GetAtomTypes().CZ,   GetAtomTypes().OH,   GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HD1,   GetAtomTypes().HD2,  GetAtomTypes().HE1,   GetAtomTypes().HE2, GetAtomTypes().HH),                                                                  GetAtomTypes().CB,  0.034, 2.940, 0.300, 6.470,  0.960,  5.660,/*Y*/   -1.0, 10.1, 10.0, 10.1, 10.0, 10.07,  0.0, 10.46, 10.00, 10.3, 10.00,  9.84, 7.50, 3.55,  9.34, 10.34,  9.11, 2.20, 0.250, 0.410,  0.112, -0.211,  0.177, -0.710,  0.700, -1.300,  0.260,  2.300, -0.210,  0.400, -0.090, -0.030, -0.038, -0.012,  0.053,  0.247, -0.015,  0.136, -0.382,  0.017, -0.114,  0.008,  0.173,  0.221, -0.445,  3.577,  1.147, -0.038, -2.688, 328.820, 6.83, 13.607,  20.23, 130.48,    1,   1, true , 1.77))), //TYR
      VAL( AddEnum( "VALINE",                /* 19 */   AATypeData( "VAL", 'V', true,  "VAL", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG1, GetAtomTypes().CG2, GetAtomTypes().H,   GetAtomTypes().HA,  GetAtomTypes().HB,   GetAtomTypes().HG11, GetAtomTypes().HG12, GetAtomTypes().HG13, GetAtomTypes().HG21, GetAtomTypes().HG22, GetAtomTypes().HG23),                                                                                                                                                                          GetAtomTypes().CB,  0.071, 3.670, 0.140, 3.000,  1.220,  6.020,/*V*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.44, 3.55,  7.44,  3.69,  9.62, 2.32, 0.270, 0.490,  0.005, -0.174,  0.243, -0.460, -2.600,  4.200,  1.080,  1.500, -1.270, -0.600, -0.170, -0.240, -0.166,  0.147,  0.066, -0.194,  0.071,  0.177, -0.044,  0.163, -0.309, -0.108,  0.311,  0.247, -2.129,  4.083,  2.297,  0.413,  5.126, 250.093, 4.27,  5.892,   0.00,  81.37,    0,   0, false, 1.87))), //VAL

      ASX( AddEnum( "ASPARAGINE_or_ASPARTIC_ACID",      AATypeData( "ASX", 'B', false, "ASP", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.052, 1.600, 0.118, 2.852, -0.698,  4.459,/*B*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,  7.46,  3.57,  9.60, 2.20, 0.233, 0.208,  0.164,  0.073, -0.181,  2.461,  7.341, -3.500, -0.849, -1.808,  0.652,  0.557,  0.289,  0.329,  0.499, -0.100, -0.185,  0.876,  0.154, -0.131,  0.156, -0.263,  0.145,  0.201, -0.205, -0.220,  0.000,  0.000,  0.000,  0.000,  0.000, 258.776, 4.16, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))), //ASX weighted average according to natural prevalence w/asn = 0.422556  asp = 0.577444
      GLX( AddEnum( "GLUTAMINE_or_GLUTAMIC_ACID",       AATypeData( "GLX", 'Z', false, "GLU", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.058, 1.560, 0.161, 3.841, -0.490,  4.006,/*Z*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,  6.96,  3.54,  9.60, 2.20, 0.399, 0.224,  0.018,  0.035, -0.046,  2.607,  6.734, -3.500, -0.779, -2.006,  0.873,  0.700,  0.189,  0.204,  0.492,  0.032, -0.276,  0.867,  0.054, -0.365,  0.100, -0.064, -0.054,  0.229,  0.128, -0.195,  0.000,  0.000,  0.000,  0.000,  0.000, 285.646, 5.22, 0.0000,  0.000,  0.000,    0,   0, false, 2.05))), //GLX weighted average according to natural prevalence w/gln = 0.35762   glu = 0.64238
      XXX( AddEnum( "ARBITRARY_AMINO_ACID",             AATypeData( "XXX", 'X', false, "XXX", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 2.165, 0.147, 3.330,  0.358,  6.141,/*X*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,  7.26,  3.57,  9.60, 2.20, 0.290, 0.298,  0.024,  0.029,  0.013,  0.624,  1.524, -0.292,  0.016, -0.042,  0.023,  0.138,  0.003,  0.008,  0.110,  0.007, -0.015,  0.241,  0.032,  0.011,  0.097,  0.047,  0.026,  0.009,  0.038,  0.020,  0.000,  0.000,  0.000,  0.000,  0.000, 263.191, 4.45, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))), //XXX weighted averaged according to natural prevalence)
      UNK( AddEnum( "UNKNOWN_AMINO_ACID",               AATypeData( "UNK", 'U', false, "UNK", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 2.165, 0.147, 3.330,  0.358,  6.141,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,  0.00,  0.00,  9.60, 2.20, 0.290, 0.298,  0.024,  0.029,  0.013,  0.624,  1.524, -0.292,  0.016, -0.042,  0.023,  0.138,  0.003,  0.008,  0.110,  0.007, -0.015,  0.241,  0.032,  0.011,  0.097,  0.047,  0.026,  0.009,  0.038,  0.020,  0.000,  0.000,  0.000,  0.000,  0.000, 263.191, 4.45, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))), //UNK weighted averaged according to natural prevalence)
      GAP( AddEnum( "SEQUENCE_GAP",                     AATypeData( "GAP", '-', false, "GAP"))),

      DAL( AddEnum( "D-ALANINE",                        AATypeData( "DAL", 'U', false, "ALA", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.085, 1.280, 0.050, 1.000,  0.310,  6.110,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.420, 0.230, -0.168,  0.097,  0.116,  0.500, -1.600,  1.800,  0.620,  0.500,  0.100, -0.300, -0.170, -0.150, -0.130,  0.075,  0.081, -0.207, -0.071, -0.136, -0.109,  0.215,  0.211, -0.044,  0.177,  0.153,  0.000,  0.000,  0.000,  0.000,  0.000, 209.020, 2.15, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      ALS( AddEnum( "2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID",
          AATypeData( "ALS", 'U', false, "ALA", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.085, 1.280, 0.050, 1.000,  0.310,  6.110,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.420, 0.230, -0.168,  0.097,  0.116,  0.500, -1.600,  1.800,  0.620,  0.500,  0.100, -0.300, -0.170, -0.150, -0.130,  0.075,  0.081, -0.207, -0.071, -0.136, -0.109,  0.215,  0.211, -0.044,  0.177,  0.153,  0.000,  0.000,  0.000,  0.000,  0.000, 209.020, 2.15, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      ACL( AddEnum( "DEOXY-CHLOROMETHYL-ARGININE",      AATypeData( "ACL", 'U', false, "ARG", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.051, 2.340, 0.290, 6.130, -1.010, 10.740,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.360, 0.250,  0.076, -0.033, -0.036,  1.810, 12.300, -4.500, -2.530, -3.000,  1.910,  1.400,  0.370,  0.320,  0.400, -0.110, -0.150,  0.978, -0.093, -0.149,  0.056, -0.084, -0.069, -0.002, -0.018, -0.068,  0.000,  0.000,  0.000,  0.000,  0.000, 335.732, 7.36, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CME( AddEnum( "S,S-(2-HYDROXYETHYL)THIOCYSTEINE", AATypeData( "CME", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CSE( AddEnum( "SELENOCYSTEINE",                   AATypeData( "CSE", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.030,  0.245, -0.148, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150,  0.040,  0.369, -0.250, -0.125,  0.240, -0.285,  1.381,  1.205, -0.159, -0.228,  0.281, -0.353,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CSD( AddEnum( "3-SULFINOALANINE",                 AATypeData( "CSD", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CSO( AddEnum( "S-HYDROXYCYSTEINE",                AATypeData( "CSO", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CSW( AddEnum( "CYSTEINE-S-DIOXIDE",               AATypeData( "CSW", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CYG( AddEnum( "2-AMINO-4-(AMINO-3-OXO-PROPYLSULFANYLCARBONYL)-BUTYRIC ACID",
          AATypeData( "CYG", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      OCS( AddEnum( "CYSTEINESULFONIC_ACID",            AATypeData( "OCS", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      SC2( AddEnum( "N-ACETYL-L-CYSTEINE",              AATypeData( "SC2", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CGU( AddEnum( "GAMMA-CARBOXY-GLUTAMIC_ACID",      AATypeData( "CGU", 'U', false, "GLU", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.069, 1.560, 0.150, 3.780, -0.640,  3.090,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.420, 0.210,  0.029,  0.058, -0.078,  3.630,  8.200, -3.500, -0.740, -3.000,  0.830,  0.700,  0.150,  0.300,  0.470,  0.100, -0.310,  0.828,  0.192, -0.352,  0.235, -0.025, -0.034, -0.184,  0.196, -0.175,  0.000,  0.000,  0.000,  0.000,  0.000, 285.025, 5.05, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      PCA( AddEnum( "PYROGLUTAMIC_ACID",                AATypeData( "PCA", 'U', false, "GLU", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.069, 1.560, 0.150, 3.780, -0.640,  3.090,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.420, 0.210,  0.029,  0.058, -0.078,  3.630,  8.200, -3.500, -0.740, -3.000,  0.830,  0.700,  0.150,  0.300,  0.470,  0.100, -0.310,  0.828,  0.192, -0.352,  0.235, -0.025, -0.034, -0.184,  0.196, -0.175,  0.000,  0.000,  0.000,  0.000,  0.000, 285.025, 5.05, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      KCX( AddEnum( "LYSINE_NZ-CARBOXYLIC ACID",        AATypeData( "KCX", 'U', false, "LYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 1.890, 0.220, 4.770, -0.990,  9.990,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.320, 0.270,  0.152, -0.043, -0.083,  2.800,  8.800, -3.900, -1.500, -3.000,  1.400,  1.800,  0.320,  0.240,  0.690, -0.130, -0.220,  1.253,  0.102, -0.129,  0.425, -0.360, -0.100, -0.025,  0.025, -0.137,  0.000,  0.000,  0.000,  0.000,  0.000, 303.428, 6.57, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      LLP( AddEnum( "2-LYSINE(3-HYDROXY-2-METHYL-5-PHOSPHONOOXYMETHYL-PYRIDIN-4-YLMETHANE)",
          AATypeData( "LLP", 'U', false, "LYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 1.890, 0.220, 4.770, -0.990,  9.990,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.320, 0.270,  0.152, -0.043, -0.083,  2.800,  8.800, -3.900, -1.500, -3.000,  1.400,  1.800,  0.320,  0.240,  0.690, -0.130, -0.220,  1.253,  0.102, -0.129,  0.425, -0.360, -0.100, -0.025,  0.025, -0.137,  0.000,  0.000,  0.000,  0.000,  0.000, 303.428, 6.57, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      M3L( AddEnum( "N-TRIMETHYLLYSINE",                AATypeData( "M3L", 'U', false, "LYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 1.890, 0.220, 4.770, -0.990,  9.990,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.320, 0.270,  0.143, -0.041, -0.078,  2.800,  8.800, -3.900, -1.500, -3.000,  1.400,  1.800,  0.320,  0.240,  0.724, -0.131, -0.220,  1.143,  0.033, -0.163,  0.402, -0.366, -0.114,  0.378, -0.085, -0.162,  0.000,  0.000,  0.000,  0.000,  0.000, 303.428, 6.57, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      MLY( AddEnum( "N-DIMETHYL-LYSINE",                AATypeData( "MLY", 'U', false, "LYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 1.890, 0.220, 4.770, -0.990,  9.990,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.320, 0.270,  0.152, -0.043, -0.083,  2.800,  8.800, -3.900, -1.500, -3.000,  1.400,  1.800,  0.320,  0.240,  0.690, -0.130, -0.220,  1.253,  0.102, -0.129,  0.425, -0.360, -0.100, -0.025,  0.025, -0.137,  0.000,  0.000,  0.000,  0.000,  0.000, 303.428, 6.57, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      LYR( AddEnum( "N~6~-[(2Z,4E,6E,8E)-3,7-DIMETHYL-9-(2,6,6-TRIMETHYLCYCLOHEX-1-EN-1-YL)NONA-2,4,6,8-TETRAENYL]LYSINE",
          AATypeData( "LYR", 'U', false, "LYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 1.890, 0.220, 4.770, -0.990,  9.990,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.320, 0.270,  0.152, -0.043, -0.083,  2.800,  8.800, -3.900, -1.500, -3.000,  1.400,  1.800,  0.320,  0.240,  0.690, -0.130, -0.220,  1.253,  0.102, -0.129,  0.425, -0.360, -0.100, -0.025,  0.025, -0.137,  0.000,  0.000,  0.000,  0.000,  0.000, 303.428, 6.57, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CXM( AddEnum( "N-CARBOXYMETHIONINE",              AATypeData( "CXM", 'U', false, "MET", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.021, 2.350, 0.220, 4.430,  1.230,  5.710,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.380, 0.320, -0.174,  0.074,  0.151, -0.670, -3.400,  1.900,  0.640,  1.300, -1.590, -0.400, -0.260, -0.190, -0.120,  0.010,  0.130, -0.255, -0.183,  0.037,  0.139,  0.441, -0.045, -0.221,  0.131,  0.354,  0.000,  0.000,  0.000,  0.000,  0.000, 291.524, 5.66, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      FME( AddEnum( "N-FORMYLMETHIONINE",               AATypeData( "FME", 'U', false, "MET", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.021, 2.350, 0.220, 4.430,  1.230,  5.710,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.380, 0.320, -0.174,  0.074,  0.151, -0.670, -3.400,  1.900,  0.640,  1.300, -1.590, -0.400, -0.260, -0.190, -0.120,  0.010,  0.130, -0.255, -0.183,  0.037,  0.139,  0.441, -0.045, -0.221,  0.131,  0.354,  0.000,  0.000,  0.000,  0.000,  0.000, 291.524, 5.66, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      MSE( AddEnum( "SELENO_METHIONINE",                AATypeData( "MSE", 'U', false, "MET", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().SE,  GetAtomTypes().CE,  GetAtomTypes().H,   GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,  GetAtomTypes().HE1,  GetAtomTypes().HE2,  GetAtomTypes().HE3),                                                                                                                                                      GetAtomTypes().CB,  0.021, 2.350, 0.220, 4.430,  1.230,  5.710,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.380, 0.320, -0.174,  0.074,  0.151, -0.670, -3.400,  1.900,  0.640,  1.300, -1.590, -0.400, -0.260, -0.190, -0.120,  0.010,  0.130, -0.255, -0.183,  0.037,  0.139,  0.441, -0.045, -0.221,  0.131,  0.354,  0.000,  0.000,  0.000,  0.000,  0.000, 291.524, 5.66, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      DPN( AddEnum( "D-PHENYLALANINE",                  AATypeData( "DPN", 'U', false, "PHE", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.039, 2.940, 0.290, 5.890,  1.790,  5.670,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.300, 0.380, -0.104, -0.016,  0.148, -1.710, -3.700,  2.800,  1.190,  2.500, -2.120, -0.500, -0.410, -0.220, -0.190, -0.020,  0.340, -0.342, -0.101,  0.371,  0.031, -0.015, -0.011,  0.097, -0.084,  0.416,  0.000,  0.000,  0.000,  0.000,  0.000, 311.302, 6.16, 0.0000,  0.000,  0.000,    0,   0, true , 2.00))),
      S1H( AddEnum( "1-HEXADECANOSULFONYL-O-L-SERINE",  AATypeData( "S1H", 'U', false, "SER", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.060, 1.310, 0.060, 1.600, -0.040,  5.700,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.200, 0.280,  0.097, -0.018, -0.068,  0.460, -0.600, -0.800, -0.180, -0.300,  0.520,  0.100,  0.050,  0.160,  0.020,  0.020, -0.040,  0.160,  0.216,  0.000, -0.032, -0.076,  0.038, -0.094, -0.097, -0.033,  0.000,  0.000,  0.000,  0.000,  0.000, 223.038, 2.66, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      SAC( AddEnum( "N-ACETYL-SERINE",                  AATypeData( "SAC", 'U', false, "SER", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.060, 1.310, 0.060, 1.600, -0.040,  5.700,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.200, 0.280,  0.097, -0.018, -0.068,  0.460, -0.600, -0.800, -0.180, -0.300,  0.520,  0.100,  0.050,  0.160,  0.020,  0.020, -0.040,  0.160,  0.216,  0.000, -0.032, -0.076,  0.038, -0.094, -0.097, -0.033,  0.000,  0.000,  0.000,  0.000,  0.000, 223.038, 2.66, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      SEP( AddEnum( "PHOSPHOSERINE",                    AATypeData( "SEP", 'U', false, "SER", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.060, 1.310, 0.060, 1.600, -0.040,  5.700,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.200, 0.280,  0.097, -0.018, -0.068,  0.460, -0.600, -0.800, -0.180, -0.300,  0.520,  0.100,  0.050,  0.160,  0.020,  0.020, -0.040,  0.160,  0.216,  0.000, -0.032, -0.076,  0.038, -0.094, -0.097, -0.033,  0.000,  0.000,  0.000,  0.000,  0.000, 223.038, 2.66, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      TPO( AddEnum( "PHOSPHOTHREONINE",                 AATypeData( "TPO", 'U', false, "THR", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.055, 3.030, 0.110, 2.600,  0.260,  5.600,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.210, 0.360,  0.104, -0.100,  0.015,  0.250, -1.200, -0.700, -0.050,  0.400,  0.070,  0.200,  0.020, -0.080, -0.010,  0.020,  0.000,  0.061,  0.136,  0.128, -0.020, -0.230, -0.092,  0.010,  0.098,  0.008,  0.000,  0.000,  0.000,  0.000,  0.000, 243.554, 3.93, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      TYS( AddEnum( "O-SULFO-L-TYROSINE",               AATypeData( "TYS", 'U', false, "TYR", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.034, 2.940, 0.300, 6.470,  0.960,  5.660,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.250, 0.410,  0.140, -0.242,  0.210, -0.710,  0.700, -1.300,  0.260,  2.300, -0.210,  0.400, -0.090, -0.030, -0.070,  0.040,  0.030,  0.259, -0.010,  0.093, -0.423,  0.005, -0.182,  0.358,  0.176,  0.141,  0.000,  0.000,  0.000,  0.000,  0.000, 328.820, 6.83, 0.0000,  0.000,  0.000,    0,   0, true , 2.00))),
      R1A( AddEnum( "methanesulfonothioate",            AATypeData( "R1A", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().SG, GetAtomTypes().SD, GetAtomTypes().CE, GetAtomTypes().C2, GetAtomTypes().C3, GetAtomTypes().C4, GetAtomTypes().C5, GetAtomTypes().C6, GetAtomTypes().C7, GetAtomTypes().C8, GetAtomTypes().C9, GetAtomTypes().N1, GetAtomTypes().O1),                                                                                                                                                      GetAtomTypes().CB,  0.000, 0.000, 0.000, 0.000,  0.000,  0.000,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.000, 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, 000.000, 0.00, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      IAS( AddEnum( "BETA-L-ASPARTIC_ACID",             AATypeData( "IAS", 'U', false, "ASP", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.058, 1.600, 0.110, 2.780, -0.770,  2.950,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.250, 0.200,  0.169,  0.067, -0.180,  3.640,  9.200, -3.500, -0.900, -3.000,  0.780,  0.600,  0.370,  0.410,  0.574, -0.068, -0.237,  0.963,  0.227, -0.193,  0.237, -0.295,  0.074,  0.266, -0.176, -0.253,  0.000,  0.000,  0.000,  0.000,  0.000, 257.993, 4.00, 0.0000,  0.000,  0.000,    0,   0, false, 2.00)))
    {
      // @see http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl
      // add mapping of atom types to pdb atom names
      ALA->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB1 , "1HB ");
      ALA->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "2HB ");
      ALA->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "3HB ");

      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG2 , "1HG ");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG3 , "2HG ");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD2 , "1HD ");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD3 , "2HD ");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HH11, "1HH1");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HH12, "2HH1");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HH21, "1HH2");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HH22, "2HH2");

      ASP->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      ASP->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      ASN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      ASN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      ASN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD21, "1HD2");
      ASN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD22, "2HD2");

      CYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      CYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      GLU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      GLU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      GLU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG2 , "1HG ");
      GLU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG3 , "2HG ");

      GLN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      GLN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      GLN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG2 , "1HG ");
      GLN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG3 , "2HG ");
      GLN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE21, "1HE2");
      GLN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE22, "2HE2");

      GLY->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HA2 , "1HA ");
      GLY->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HA3 , "2HA ");

      HIS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      HIS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG12, "1HG1");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG13, "2HG1");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG21, "1HG2");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG22, "2HG2");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG23, "3HG2");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD11, "1HD1");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD12, "2HD1");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD13, "3HD1");

      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD11, "1HD1");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD12, "2HD1");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD13, "3HD1");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD21, "1HD2");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD22, "2HD2");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD23, "3HD2");

      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG2 , "1HG ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG3 , "2HG ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD2 , "1HD ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD3 , "2HD ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE2 , "1HE ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE3 , "2HE ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HZ1 , "1HZ ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HZ2 , "2HZ ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HZ3 , "3HZ ");

      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG2 , "1HG ");
      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG3 , "2HG ");
      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE1 , "1HE ");
      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE2 , "2HE ");
      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE3 , "3HE ");

      PHE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      PHE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      PRO->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      PRO->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      PRO->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG2 , "1HG ");
      PRO->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG3 , "2HG ");
      PRO->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD2 , "1HD ");
      PRO->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD3 , "2HD ");

      SER->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      SER->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      THR->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG21, "1HG2");
      THR->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG22, "2HG2");
      THR->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG23, "3HG2");

      TRP->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      TRP->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      TYR->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      TYR->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      VAL->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG11, "1HG1");
      VAL->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG12, "2HG1");
      VAL->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG13, "3HG1");
      VAL->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG21, "1HG2");
      VAL->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG22, "2HG2");
      VAL->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG23, "3HG2");

      // set pka for commonly charged hydrogens
      ARG->SetSideChainIonizableHType( GetAtomTypes().HD22);
      HIS->SetSideChainIonizableHType( GetAtomTypes().HD1);
      LYS->SetSideChainIonizableHType( GetAtomTypes().HZ3);
      ASP->SetSideChainIonizableHType( GetAtomTypes().HD2);
      GLU->SetSideChainIonizableHType( GetAtomTypes().HE2);
      CYS->SetSideChainIonizableHType( GetAtomTypes().HG);
      TYR->SetSideChainIonizableHType( GetAtomTypes().HH);

      // set the dihedral angle atoms
      MET->SetDihedralAngles
      (
        storage::Map< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >::Create
        (
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_One,
            storage::VectorND< 4, AtomType>( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().CB, GetAtomTypes().CG)
          ),
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_Two,
            storage::VectorND< 4, AtomType>( GetAtomTypes().CA, GetAtomTypes().CB, GetAtomTypes().CG, GetAtomTypes().SD)
          ),
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_Three,
            storage::VectorND< 4, AtomType>( GetAtomTypes().CB, GetAtomTypes().CG, GetAtomTypes().SD, GetAtomTypes().CE)
          )
        )
      );

      R1A->SetDihedralAngles
      (
        storage::Map< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >::Create
        (
          std::make_pair
          (
            ChiAngle::e_One,
            storage::VectorND< 4, AtomType>( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().CB, GetAtomTypes().SG)
          ),
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_Two,
            storage::VectorND< 4, AtomType>( GetAtomTypes().CA, GetAtomTypes().CB, GetAtomTypes().SG, GetAtomTypes().SD)
          ),
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_Three,
            storage::VectorND< 4, AtomType>( GetAtomTypes().CB, GetAtomTypes().SG, GetAtomTypes().SD, GetAtomTypes().CE)
          ),
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_Four,
            storage::VectorND< 4, AtomType>( GetAtomTypes().SG, GetAtomTypes().SD, GetAtomTypes().CE, GetAtomTypes().C3)
          ),
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_Five,
            storage::VectorND< 4, AtomType>( GetAtomTypes().SD, GetAtomTypes().CE, GetAtomTypes().C3, GetAtomTypes().C4)
          )
        )
      );

      // determine the parent types
      for( iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        ( *itr)->DetermineParentType( *this);
      }
      e_Undefined->DetermineParentType( *this);
    }

    //! @brief function to deduce AAType from three letter code of an amino acid
    //! @param THREE_LETTER_CODE three letter code descriptor for amino acid of interest
    //! @return AAType specified by given THREE_LETTER_CODE
    const AAType &AATypes::AATypeFromThreeLetterCode( const std::string &THREE_LETTER_CODE) const
    {
      // assert the provided string has correct size
      if( THREE_LETTER_CODE.size() != 3)
      {
        BCL_MessageVrb( "three letter code should consist of three letters: " + THREE_LETTER_CODE);
        return e_Undefined;
      }

      // iterate over all AATypes
      for
      (
        const_iterator type_itr( GetAATypes().Begin()), type_itr_end( AATypes::End());
        type_itr != type_itr_end;
        ++type_itr
      )
      {
        // if THREE_LETTER_CODE matches this AATypes's three letter code
        if( ( *type_itr)->GetThreeLetterCode() == THREE_LETTER_CODE)
        {
          // return this AAType
          return *type_itr;
        }
      }

      // if no match was found return undefined AAType
      return GetAATypes().e_Undefined;
    }

    //! @brief return if two amino acids have same parent
    //! this function returns true if two amino acids denoted by there three letter code have the same parent.
    //! MSE == MET would return true
    //! @param THREE_LETTER_CODE_LHS three letter code descriptor for amino acid of interest lhs
    //! @param THREE_LETTER_CODE_RHS three letter code descriptor for amino acid of interest rhs
    //! @return true if the two amino acids have same parent
    bool AATypes::HaveSameParent( const std::string &THREE_LETTER_CODE_LHS, const std::string &THREE_LETTER_CODE_RHS) const
    {
      if( THREE_LETTER_CODE_LHS == THREE_LETTER_CODE_RHS)
      {
        return true;
      }

      const AAType &lhs( AATypeFromThreeLetterCode( THREE_LETTER_CODE_LHS));
      const AAType &rhs( AATypeFromThreeLetterCode( THREE_LETTER_CODE_RHS));
      if( !lhs.IsDefined() || !rhs.IsDefined())
      {
        return false;
      }
      return lhs->GetParentType() == rhs->GetParentType();
    }

    namespace
    {
      // anonymous namespace to prevent unwanted symbol export
      //! @brief create a vector of 256 values that return, for a given character, the corresponding aa
      //! @return a vector of 256 values that return, for a given character, the corresponding aa
      storage::Vector< AAType> CreateAATypeFromOneLetterCodeVector()
      {
        storage::Vector< AAType> aatypes( 256, GetAATypes().e_Undefined);
        // iterate over all AATypes
        for
        (
          AATypes::const_iterator type_itr( GetAATypes().Begin()), type_itr_end( GetAATypes().End());
          type_itr != type_itr_end;
          ++type_itr
        )
        {
          // check that the aatypes vector was not already set; typically this happens for unnatural types,
          // usually denoted with U as the one letter code, but we want to match
          if( !aatypes( size_t( ( *type_itr)->GetOneLetterCode())).IsDefined())
          {
            aatypes( size_t( ( *type_itr)->GetOneLetterCode())) = *type_itr;
          }
        }
        return aatypes;
      }
    }

    //! @brief function to deduce AAType from one letter code of an amino acid
    //! @param ONE_LETTER_CODE one letter code descriptor for amino acid of interest
    //! @return AAType specified by the given ONE_LETTER_CODE
    const AAType &AATypes::AATypeFromOneLetterCode( const char ONE_LETTER_CODE) const
    {
      static const storage::Vector< AAType> s_aa_type_one_letter_code( CreateAATypeFromOneLetterCodeVector());
      // return the type from the vector
      return size_t( ONE_LETTER_CODE) < s_aa_type_one_letter_code.GetSize()
             ? s_aa_type_one_letter_code( size_t( ONE_LETTER_CODE))
             : GetAATypes().e_Undefined;
    }

    //! @brief gives the 20 natural amino acid types
    //! @return set of the 20 natural amino acid types
    const storage::Set< AAType> &AATypes::GetNaturalAATypes() const
    {
      static storage::Set< AAType> s_natural_aa_types;

      if( s_natural_aa_types.IsEmpty())
      {
        // iterate over AATypes
        for
        (
          const_iterator type_itr( Begin()), type_itr_end( ASX.GetIterator());
          type_itr != type_itr_end;
          ++type_itr
        )
        {
          s_natural_aa_types.Insert( *type_itr);
        }
        BCL_Assert
        (
          20 == s_natural_aa_types.GetSize(),
          "number natural residue types should be " + util::Format()( 20)
          + " but is " + util::Format()( s_natural_aa_types.GetSize())
        );
      }

      return s_natural_aa_types;
    }

    //! @brief construct on access function for all AATypes
    //! @return reference to only instance of AATypes enum
    const AATypes &GetAATypes()
    {
      return AATypes::GetEnums();
    }

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< biol::AATypeData, biol::AATypes>;

  } // namespace util
} // namespace bcl
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
#include "biol/bcl_biol_aa_compare.h"
#include "biol/bcl_biol_align_by_aa_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AlignByAAData::s_Instance
    (
      GetObjectInstances().AddInstance( new AlignByAAData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AlignByAAData::AlignByAAData()
    {
    }

    //! @brief Clone function
    //! @return pointer to new AlignByAAData
    AlignByAAData *AlignByAAData::Clone() const
    {
      return new AlignByAAData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AlignByAAData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Align is the function which does the actually aligns data together
    //! @param ALIGNMENT_A is the Alignment to be aligned with ALIGNMENT_B
    //! @param ALIGNMENT_B is the Alignment to be aligned with ALIGNMENT_A
    //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
    storage::Pair< align::AlignmentNode< AABase>, double> AlignByAAData::AlignPair
    (
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_A,
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_B
    ) const
    {
      // end
      return AlignPair( ALIGNMENT_A, ALIGNMENT_B, AACompareDataPtr());
    }

    //! @brief Align is the function which does the actually aligns data together
    //! @param ALIGNMENT_A is the Alignment to be aligned with ALIGNMENT_B
    //! @param ALIGNMENT_B is the Alignment to be aligned with ALIGNMENT_A
    //! @param COMPARISON method for comparing the aa data to see if it is equal or not
    //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
    storage::Pair< align::AlignmentNode< AABase>, double> AlignByAAData::AlignPair
    (
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_A,
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_B,
      const util::BinaryFunctionInterface< AABase, AABase, bool> &COMPARISON
    ) const
    {
      // initialize return type
      align::AlignmentNode< AABase> node( ALIGNMENT_A, ALIGNMENT_B);
      const double score( AlignPairWithNode( node, COMPARISON));
      return storage::Pair< align::AlignmentNode< AABase>, double>( node, score);
    }

    //! @brief Align is the function which does the actually aligns data together
    //! @param NODE is the alignment node, which must alread contain the two parent alignment interfaces
    //! @param COMPARISON method for comparing the aa data to see if it is equal or not
    //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
    double AlignByAAData::AlignPairWithNode
    (
      align::AlignmentInterface< AABase> &NODE,
      const util::BinaryFunctionInterface< AABase, AABase, bool> &COMPARISON
    ) const
    {
      iterate::Generic< const align::AlignmentInterface< AABase> >
        itr_child_alignments( NODE.GetChildAlignmentsIterator());
      BCL_Assert
      (
        itr_child_alignments.GetSize() == size_t( 2),
        "Expected two child alignments already in the node"
      );
      const align::AlignmentInterface< AABase> &alignment_first( *itr_child_alignments);
      const align::AlignmentInterface< AABase> &alignment_second( *++itr_child_alignments);

      const align::AlignmentInterface< AABase> &alignment_a
      (
        alignment_first.GetSize() >= alignment_second.GetSize()
        ? alignment_first
        : alignment_second
      );
      const align::AlignmentInterface< AABase> &alignment_b
      (
        alignment_first.GetSize() >= alignment_second.GetSize()
        ? alignment_second
        : alignment_first
      );

      // make sure both have depth of at least 1
      if( alignment_a.GetDepth() == 0 || alignment_b.GetDepth() == 0)
      {
        BCL_MessageCrt( "The provided alignments should have depths of at least 1!");
        return util::GetUndefinedDouble();
      }
      double score( util::GetUndefinedDouble());

      // initialize iterators
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_a( alignment_a.GetAssignments().Begin());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_a_end( alignment_a.GetAssignments().End());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_b( alignment_b.GetAssignments().Begin());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_b_end( alignment_b.GetAssignments().End());

      // initialize empty pointer to be used for gaps
      const util::SiPtr< const AABase> empty_ptr;

      // while there are still elements left in ALIGNMENT_B
      while( itr_b != itr_b_end)
      {
        // check that a is not at the end yet
        if( itr_a == itr_a_end)
        {
          BCL_MessageCrt( "The provided alignment a should have all aadata of b!");
          NODE.ResetAssignments();
          score = util::GetUndefined< double>();
          return score;
        }

        const align::Assignment< AABase> &assignment_a( **itr_a), &assignment_b( **itr_b);
        BCL_MessageVrb( "Trying to align itr_a=" + assignment_a.ToString() + " and itr_b=" + assignment_b.ToString());

        // iterate until the AAData pointer by both pointers are the same
        if( !COMPARISON( *( assignment_a.GetMembers().FirstElement()), *( assignment_b.GetMembers().FirstElement())))
        {
          BCL_MessageDbg( "Assignment comparison returned false, not aligning");

          // insert assignment of itr_a with gap
          NODE.Append
          (
            util::ShPtr< align::Assignment< AABase> >
            (
              new align::Assignment< AABase>( assignment_a.GetMembers().FirstElement(), empty_ptr)
            )
          );

          // move to next amino acid
          ++itr_a;
        }
        else // aa data pointer agrees
        {
          BCL_MessageDbg( "Assignment comparison returned true, aligning");

          // insert assignment of itr_a and itr_b
          NODE.Append
          (
            util::ShPtr< align::Assignment< AABase> >
            (
              new align::Assignment< AABase>
              (
                assignment_a.GetMembers().FirstElement(), assignment_b.GetMembers().FirstElement()
              )
            )
          );

          // move both iterators
          ++itr_a;
          ++itr_b;
        }
      } // end of while loop

      // iterate while itr_a reaches end
      while( itr_a != itr_a_end)
      {
        // insert assignment of itr_a with gap
        NODE.Append
        (
          util::ShPtr< align::Assignment< AABase> >
          (
            new align::Assignment< AABase>( ( *itr_a)->GetMembers().FirstElement(), empty_ptr)
          )
        );

        // move the next residue
        ++itr_a;
      }

      // end
      return score;
    }

    //! @brief Align is the function which does the actually aligns data together
    //! @param ALIGNMENT_A is the Alignment to be aligned with ALIGNMENT_B
    //! @param ALIGNMENT_B is the Alignment to be aligned with ALIGNMENT_A
    //! @param TEMPLATE_ALIGNMENT is the template alignment, ALIGNMENT_A and ALIGNMENT_B will be aligned according to this
    //!        alignment - the first sequence must correspond to ALIGNMENT_A and the second to ALIGNMENT_B
    //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
    storage::Pair< align::AlignmentNode< AABase>, double> AlignByAAData::AlignPair
    (
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_A,
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_B,
      const util::ShPtr< align::AlignmentInterface< AABase> > &TEMPLATE_ALIGNMENT
    ) const
    {
      // initialize return type
      align::AlignmentNode< AABase> alignment( ALIGNMENT_A, ALIGNMENT_B);
      const double score( AlignPairWithNode( alignment, TEMPLATE_ALIGNMENT));
      return storage::Pair< align::AlignmentNode< AABase>, double>( alignment, score);
    }

    //! @brief Align is the function which does the actually aligns data together
    //! @param NODE is the alignment node, which must alread contain the two parent alignment interfaces
    //! @param COMPARISON method for comparing the aa data to see if it is equal or not
    //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
    double AlignByAAData::AlignPairWithNode
    (
      align::AlignmentInterface< AABase> &NODE,
      const util::ShPtr< align::AlignmentInterface< AABase> > &TEMPLATE_ALIGNMENT
    ) const
    {
      iterate::Generic< const align::AlignmentInterface< AABase> >
        itr_child_alignments( NODE.GetChildAlignmentsIterator());
      BCL_Assert
      (
        itr_child_alignments.GetSize() == size_t( 2),
        "Expected two child alignments already in the node"
      );
      const align::AlignmentInterface< AABase> &alignment_a( *itr_child_alignments);
      const align::AlignmentInterface< AABase> &alignment_b( *++itr_child_alignments);

      double score( util::GetUndefinedDouble());

      // initialize iterators
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_a( alignment_a.GetAssignments().Begin());
      const util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_a_end( alignment_a.GetAssignments().End());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_b( alignment_b.GetAssignments().Begin());
      const util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_b_end( alignment_b.GetAssignments().End());

      // initialize empty pointer to be used for gaps
      const util::SiPtr< const AABase> empty_ptr;

      // iterate through the template alignment
      for
      (
        util::ShPtrList< align::Assignment< AABase> >::const_iterator
          template_itr( TEMPLATE_ALIGNMENT->GetAssignments().Begin()),
          template_itr_end( TEMPLATE_ALIGNMENT->GetAssignments().End());
        template_itr != template_itr_end; ++template_itr
      )
      {
        // get the AA's
        util::SiPtr< const AABase> sp_aa_a( empty_ptr);
        util::SiPtr< const AABase> sp_aa_b( empty_ptr);
        const util::SiPtr< const AABase> &sp_aa_template_a( ( *template_itr)->GetMembers().FirstElement());
        const util::SiPtr< const AABase> &sp_aa_template_b( ( *template_itr)->GetMembers().LastElement());

        // if template is >= alignment A
        if
        (
          itr_a != itr_a_end &&
          sp_aa_template_a.IsDefined() &&
          !AALessThanSeqID()( sp_aa_template_a, ( *itr_a)->GetMembers().FirstElement())
        )
        {
          // iterate along the alignment until it is no longer smaller
          while( itr_a != itr_a_end && AALessThanSeqID()( ( *itr_a)->GetMembers().FirstElement(), sp_aa_template_a))
          {
            ++itr_a;
          }

          // if the data match
          if( itr_a != itr_a_end && AACompareData()( sp_aa_template_a, ( *itr_a)->GetMembers().FirstElement()))
          {
            // use this assignment
            sp_aa_a = ( *itr_a)->GetMembers().FirstElement();

            // increment
            ++itr_a;
          }
        }

        // if template is >= alignment B
        if
        (
          itr_b != itr_b_end &&
          sp_aa_template_b.IsDefined() &&
          !AALessThanSeqID()( sp_aa_template_b, ( *itr_b)->GetMembers().FirstElement())
        )
        {
          // iterate along the alignment until it is no longer smaller
          while( itr_b != itr_b_end && AALessThanSeqID()( ( *itr_b)->GetMembers().FirstElement(), sp_aa_template_b))
          {
            ++itr_b;
          }

          // if the data match
          if( itr_b != itr_b_end && AACompareData()( sp_aa_template_b, ( *itr_b)->GetMembers().FirstElement()))
          {
            // use this assignment
            sp_aa_b = ( *itr_b)->GetMembers().FirstElement();

            // increment
            ++itr_b;
          }
        }

        // pushback the assignment
        NODE.Append
        (
          util::ShPtr< align::Assignment< AABase> >
          (
            new align::Assignment< AABase>( sp_aa_a, sp_aa_b)
          )
        );
      }
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AlignByAAData::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AlignByAAData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_align_by_pdb_id.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AlignByPdbID::s_Instance
    (
      GetObjectInstances().AddInstance( new AlignByPdbID())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AlignByPdbID::AlignByPdbID()
    {
    }

    //! @brief Clone function
    //! @return pointer to new AlignByPdbID
    AlignByPdbID *AlignByPdbID::Clone() const
    {
      return new AlignByPdbID( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AlignByPdbID::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Align is the function which does the actually aligns data together
    //! @param ALIGNMENT_A is the Alignment to be aligned with ALIGNMENT_B
    //! @param ALIGNMENT_B is the Alignment to be aligned with ALIGNMENT_A
    //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
    storage::Pair< align::AlignmentNode< AABase>, double> AlignByPdbID::AlignPair
    (
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_A,
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_B
    ) const
    {
      // make sure both have depth of at least 1
      if( ALIGNMENT_A->GetDepth() == 0 || ALIGNMENT_B->GetDepth() == 0)
      {
        BCL_MessageCrt( "The provided alignments should have depths of at least 1!");
        return storage::Pair< align::AlignmentNode< AABase>, double>( align::AlignmentNode< AABase>(), util::GetUndefinedDouble());
      }

      // initialize return type
      align::AlignmentNode< AABase> alignment( ALIGNMENT_A, ALIGNMENT_B);
      storage::Pair< align::AlignmentNode< AABase>, double> alignment_and_score( alignment, util::GetUndefinedDouble());

      // initialize iterators
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_a( ALIGNMENT_A->GetAssignments().Begin());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_a_end( ALIGNMENT_A->GetAssignments().End());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_b( ALIGNMENT_B->GetAssignments().Begin());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_b_end( ALIGNMENT_B->GetAssignments().End());

      // initialize empty pointer to be used for gaps
      const util::SiPtr< const AABase> empty_ptr;

      // while there are still elements left in ALIGNMENT_B
      while( itr_a != itr_a_end && itr_b != itr_b_end)
      {
        util::SiPtr< const AABase> amino_acid_a( ( *itr_a)->GetMembers().FirstElement());
        util::SiPtr< const AABase> amino_acid_b( ( *itr_b)->GetMembers().FirstElement());

        // match the chains
        if( amino_acid_a->GetChainID() == amino_acid_b->GetChainID())
        {
          // match pdb id
          if( amino_acid_a->GetPdbID() == amino_acid_b->GetPdbID())
          {
            // match pdb icode
            if( amino_acid_a->GetPdbICode() == amino_acid_b->GetPdbICode())
            {
              // move both iterators
              ++itr_a;
              ++itr_b;
            }
            // pdb i code smaller
            else if( amino_acid_a->GetPdbICode() < amino_acid_b->GetPdbICode())
            {
              ++itr_a;
              amino_acid_b = empty_ptr;
            }
            // pdb i code larger
            else
            {
              ++itr_b;
              amino_acid_a = empty_ptr;
            }
          }
          // pdb id smaller
          else if( amino_acid_a->GetPdbID() < amino_acid_b->GetPdbID())
          {
            ++itr_a;
            amino_acid_b = empty_ptr;
          }
          // pdb id larger
          else
          {
            ++itr_b;
            amino_acid_a = empty_ptr;
          }
        }
        // chain id smaller
        else if( amino_acid_a->GetChainID() < amino_acid_b->GetChainID())
        {
          ++itr_a;
          amino_acid_b = empty_ptr;
        }
        // chain id larger
        else
        {
          ++itr_b;
          amino_acid_a = empty_ptr;
        }

        // insert assignment of itr_a and itr_b
        alignment_and_score.First().Append
        (
          util::ShPtr< align::Assignment< AABase> >( new align::Assignment< AABase>( amino_acid_a, amino_acid_b))
        );
      } // end of while loop

      // iterate while itr_a reaches end
      while( itr_a != itr_a_end)
      {
        // insert assignment of itr_a with gap
        alignment_and_score.First().Append
        (
          util::ShPtr< align::Assignment< AABase> >
          (
            new align::Assignment< AABase>( ( *itr_a)->GetMembers().FirstElement(), empty_ptr)
          )
        );

        // move the next residue
        ++itr_a;
      }

      // iterate while itr_b reaches end
      while( itr_b != itr_b_end)
      {
        // insert assignment of itr_b with gap
        alignment_and_score.First().Append
        (
          util::ShPtr< align::Assignment< AABase> >
          (
            new align::Assignment< AABase>( empty_ptr, ( *itr_b)->GetMembers().FirstElement())
          )
        );

        // move the next residue
        ++itr_b;
      }

      // end
      return alignment_and_score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AlignByPdbID::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AlignByPdbID::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_atom.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Atom::s_Instance
    (
      GetObjectInstances().AddInstance( new Atom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct undefined atom
    Atom::Atom() :
      m_Coordinates( util::GetUndefined< double>()),
      m_Type( util::GetUndefined< AtomType>()),
      m_PdbID( util::GetUndefined< size_t>()),
      m_BFactor( util::GetUndefined< double>())
    {
    }

    //! @brief construct atom from a single given ATOM_TYPE
    //! @param ATOM_TYPE AtomType of interest
    Atom::Atom( const AtomType &ATOM_TYPE) :
      m_Coordinates( util::GetUndefined< double>()),
      m_Type( ATOM_TYPE),
      m_PdbID( util::GetUndefined< size_t>()),
      m_BFactor( util::GetUndefined< double>())
    {
    }

    //! @brief construct atom from atom type ( string name can also be used instead), pdb id and b factor
    //! @param ATOM_TYPE AtomType of interest
    //! @param PDB_ID pdb id of atom
    //! @param B_FACTOR B factor of atom
    Atom::Atom
    (
      const AtomType &ATOM_TYPE,
      const size_t PDB_ID,
      const double B_FACTOR
    ) :
      m_Coordinates( util::GetUndefined< double>()),
      m_Type( ATOM_TYPE),
      m_PdbID( PDB_ID),
      m_BFactor( B_FACTOR)
    {
    }

    //! @brief construct atom from coordinates, atom type (string name can also be used instead), pdb id and b factor
    //! @param COORDINATES 3D coordinates of the atom
    //! @param ATOM_TYPE AtomType of interest
    //! @param PDB_ID pdb id of atom
    //! @param B_FACTOR B factor of atom
    Atom::Atom
    (
      const linal::Vector3D &COORDINATES,
      const AtomType &ATOM_TYPE,
      const size_t PDB_ID,
      const double B_FACTOR
    ) :
      m_Coordinates( COORDINATES),
      m_Type( ATOM_TYPE),
      m_PdbID( PDB_ID),
      m_BFactor( B_FACTOR)
    {
    }

    //! @brief virtual copy constructor
    Atom *Atom::Clone() const
    {
      return new Atom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Atom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! return coordinates as Vector3D
    const linal::Vector3D &Atom::GetCoordinates() const
    {
      return m_Coordinates;
    }

    //! change coordinates and increase state
    void Atom::SetCoordinates( const linal::Vector3D &COORDINATES)
    {
      m_Coordinates = COORDINATES;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief static function to find an atom of the specified ATOM_TYPE from given ATOM_LIST
    //! @param ATOM_LIST SiPtrVector of Atoms
    //! @param ATOM_TYPE AtomType of interest
    //! @return SiPtr to requested ATOM_TYPE from given ATOM_LIST  return undefined SiPtr if not found
    const util::SiPtr< const Atom> &Atom::FindAtom
    (
      const util::SiPtrVector< const Atom> &ATOM_LIST,
      const AtomType &ATOM_TYPE
    )
    {
      // initialize undefined atom
      static const Atom s_undefined_atom;

      // initialize pointer to undefined atom
      static const util::SiPtr< const Atom> s_undefined_atom_pointer( s_undefined_atom);

      // iterate over atoms in the list
      for
      (
        util::SiPtrVector< const Atom>::const_iterator atom_itr( ATOM_LIST.Begin()), atom_itr_end( ATOM_LIST.End());
        atom_itr != atom_itr_end;
        ++atom_itr
      )
      {
        // if the type of this atom matches the specified type, return it
        if( ( *atom_itr)->GetType() == ATOM_TYPE)
        {
          return *atom_itr;
        }
      }

      // if the atom was not found return the simple pointer to undefined atom
      return s_undefined_atom_pointer;
    }

    //! transforms the coordinates with a transformationmatrix3D
    void Atom::Transform( const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D)
    {
      m_Coordinates.Transform( TRANSFORMATIONMATRIX3D);
    }

    //! translate coordinates
    void Atom::Translate( const linal::Vector3D &TRANSLATION)
    {
      m_Coordinates.Translate( TRANSLATION);
    }

    //! rotates coordinates
    void Atom::Rotate( const math::RotationMatrix3D &ROTATIONMATRIX3D)
    {
      m_Coordinates.Rotate( ROTATIONMATRIX3D);
    }

    //! returns the Center
    linal::Vector3D Atom::GetCenter() const
    {
      return m_Coordinates;
    }

    //! @brief AllCoordinatesDefined determines if the x,y,and z coordinates for the atom are defined
    //! @return returns a bool - true if all coordinates are defined, false otherwise
    bool Atom::AllCoordinatesDefined() const
    {
      return m_Coordinates.IsDefined();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read Atom from io::IFStream
    std::istream &Atom::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Coordinates, ISTREAM);
      io::Serialize::Read( m_Type       , ISTREAM);
      io::Serialize::Read( m_PdbID      , ISTREAM);
      io::Serialize::Read( m_BFactor    , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &Atom::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Coordinates, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Type       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PdbID      , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BFactor    , OSTREAM, INDENT);

      // end and return the ostream
      return OSTREAM;
    }

  ////////////////
  // operations //
  ////////////////

    //! returns square distance between two atoms
    double SquareDistance( const Atom &ATOM_A, const Atom &ATOM_B)
    {
      return linal::SquareDistance
      (
        ATOM_A.GetCoordinates(),
        ATOM_B.GetCoordinates()
      );
    }

    //! returns distance between two atoms
    double Distance( const Atom &ATOM_A, const Atom &ATOM_B)
    {
      return linal::Distance
      (
        ATOM_A.GetCoordinates(),
        ATOM_B.GetCoordinates()
      );
    }

    //! returns angle between three atoms
    double ProjAngle( const Atom &ATOM_A, const Atom &ATOM_B, const Atom &ATOM_C)
    {
      return linal::ProjAngle
      (
        ATOM_A.GetCoordinates(),
        ATOM_B.GetCoordinates(),
        ATOM_C.GetCoordinates()
      );
    }

    //! returns dihedral between four atoms
    double Dihedral( const Atom &ATOM_A, const Atom &ATOM_B, const Atom &ATOM_C, const Atom &ATOM_D)
    {
      return linal::Dihedral
      (
        ATOM_A.GetCoordinates(),
        ATOM_B.GetCoordinates(),
        ATOM_C.GetCoordinates(),
        ATOM_D.GetCoordinates()
      );
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_atom_group_type_data.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom_types.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_sum_function.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AtomGroupTypeData::s_Instance
    (
      GetObjectInstances().AddInstance( new AtomGroupTypeData())
    );

    //! solvent density
    const double AtomGroupTypeData::s_SolventDensity( 0.334);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct undefined atom type
    AtomGroupTypeData::AtomGroupTypeData() :
      m_Vacuo(),
      m_Mass( util::GetUndefined< double>()),
      m_SurfaceArea( util::GetUndefined< double>())
    {
    }

    //! @brief construct atom type from given information
    //! @param BASE_ATOM_TYPE base atom type
    //! @param DISPLACED_SOLVENT_VOLUME Volume of solvent displaced by atom group
    //! @param RADIUS Radius of theoretical sphere of displaced solvent
    //! @param HYDROGEN_COUNT Number of Hydrogens bound to Heavy Atom
    AtomGroupTypeData::AtomGroupTypeData
    (
      const AtomType &BASE_ATOM_TYPE,
      const double &DISPLACED_SOLVENT_VOLUME,
      const double &RADIUS,
      const size_t &HYDROGEN_COUNT,
      const double &H2O_SCATTERING_LENGTH,
      const double &D2O_SCATTERING_LENGTH
    ) :
      m_Vacuo(),
      m_Water(),
      m_BaseAtomType( BASE_ATOM_TYPE.GetName()),
      m_DisplacedSolventVolume( DISPLACED_SOLVENT_VOLUME),
      m_Radius( RADIUS),
      m_HydrogenCount( HYDROGEN_COUNT),
      m_H2OScatteringLength( H2O_SCATTERING_LENGTH),
      m_D2OScatteringLength( D2O_SCATTERING_LENGTH),
      m_Mass( DISPLACED_SOLVENT_VOLUME * s_SolventDensity),
      m_SurfaceArea( ( -math::Pow( DISPLACED_SOLVENT_VOLUME, 2.0 / 3.0)) / ( 4.0 * math::g_Pi))
    {

      // initialize structure factors for water
      util::ShPtr< math::SumFunction< restraint::SasDataParameters, double> > water_factors
      (
        new math::SumFunction< restraint::SasDataParameters, double>()
      );

      // get the crommer mann constants for hydrogen
      static const math::SumFunction< restraint::SasDataParameters, double>
      h_form_factor( GetAtomTypes().H->GetElementType()->GetStructureFactor(), 1.0, 0.0);

      // get the crommer mann constants for oxygen
      static const math::SumFunction< restraint::SasDataParameters, double>
      o_form_factor( GetAtomTypes().O->GetElementType()->GetStructureFactor(), 1.0, 0.0);

      // add hydrogen structure factors to water_factors
      *water_factors += double( 2.0) * h_form_factor;

      // add oxygen structure factors to water_factors
      *water_factors += double( 1.0) * o_form_factor;
      m_Water = *water_factors;

      // initialize structure factors
      math::SumFunction< restraint::SasDataParameters, double> structure_factors;

      // initialize with element structure factor
      structure_factors += BASE_ATOM_TYPE->GetElementType()->GetStructureFactor();

      // add hydrogen structure factors (use h_form_factor from above)
      structure_factors += double( HYDROGEN_COUNT) * h_form_factor;

      // set vacuo form factor
      m_Vacuo = structure_factors;
    }

    //! @brief virtual copy constructor
    AtomGroupTypeData *AtomGroupTypeData::Clone() const
    {
      return new AtomGroupTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomGroupTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief () operator
    //! @param FORM_FACTOR_DATA Data to Calculate the form factors
    //! @return form factor
    double AtomGroupTypeData::operator()( const restraint::SasDataParameters &FORM_FACTOR_DATA) const
    {

      const double &Qvalue( FORM_FACTOR_DATA.GetQValue());
      const double &Sasa( FORM_FACTOR_DATA.GetSasaValue() / 100);
      const double &ExcludedVolumeParameter( FORM_FACTOR_DATA.GetExcludedVolume());
      const double &HydrationShellParameter( FORM_FACTOR_DATA.GetHydrationShell());
      const double &DeteriumPercentage( FORM_FACTOR_DATA.GetDeuteriumExchangeRate());

      double FormFactor( util::GetUndefinedDouble());

      if( FORM_FACTOR_DATA.GetSansImplementation())
      {
        double water_contribution( -0.00562);
        double deuterium_contribution( 0.0697 * FORM_FACTOR_DATA.GetDeuteriumExchangeRate());
        double solvent_scattering_length_density( water_contribution + deuterium_contribution);

        double vacuumScatteringLength( 0.0);

        // Hydrodens bound to Carbon will not exchange
        if( m_BaseAtomType == "C")
        {
          vacuumScatteringLength = m_H2OScatteringLength;
        }

        // Perform the deuterium exchange for Nitrogen, Oxygen, and Sulfer
        if( m_BaseAtomType == "N" || m_BaseAtomType == "O" || m_BaseAtomType == "S")
        {
          vacuumScatteringLength = ( m_H2OScatteringLength * ( 1 - DeteriumPercentage)) + ( m_D2OScatteringLength * DeteriumPercentage);
        }

        if( m_HydrogenCount == 0)
        {
          vacuumScatteringLength = m_H2OScatteringLength;
        }

        // Perform Calculation of water factor based on Deterium Exchange rate

        double water( -0.168);
        double deuterium( 1.915);

        double hydration_scattering_length( water * ( 1 - DeteriumPercentage) + ( deuterium * DeteriumPercentage));

        FormFactor =
          vacuumScatteringLength -
          ExcludedVolumeParameter *solvent_scattering_length_density * m_DisplacedSolventVolume * std::exp( math::Sqr( Qvalue) * m_SurfaceArea) +
          HydrationShellParameter *Sasa *hydration_scattering_length;
      }
      else
      {
        FormFactor =
          m_Vacuo->operator ()( FORM_FACTOR_DATA) -
          ExcludedVolumeParameter *m_Mass * std::exp( math::Sqr( Qvalue) * m_SurfaceArea) +
          HydrationShellParameter *Sasa *m_Water->operator()( FORM_FACTOR_DATA);
      }

      return FormFactor;
    }
    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomGroupTypeData::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Stores biological and group properties of atom groups.");
      serializer.AddInitializer
      (
        "base atom type",
        "",
        io::Serialization::GetAgent( &m_BaseAtomType)
      );
      serializer.AddInitializer
      (
        "displaced solvent volume",
        "",
        io::Serialization::GetAgent( &m_DisplacedSolventVolume)
      );
      serializer.AddInitializer
      (
        "radius",
        "",
        io::Serialization::GetAgent( &m_Radius)
      );
      serializer.AddInitializer
      (
        "hydrogen count",
        "",
        io::Serialization::GetAgent( &m_HydrogenCount)
      );
      serializer.AddInitializer
      (
        "h2o scattering length",
        "",
        io::Serialization::GetAgent( &m_H2OScatteringLength)
      );
      serializer.AddInitializer
      (
        "s2o scattering length",
        "",
        io::Serialization::GetAgent( &m_D2OScatteringLength)
      );
      serializer.AddDataMember
      (
        "mass",
        io::Serialization::GetAgent( &m_Mass)
      );
      serializer.AddDataMember
      (
        "surface area",
        io::Serialization::GetAgent( &m_SurfaceArea)
      );
      serializer.AddDataMember
      (
        "form factor vacuo",
        io::Serialization::GetAgent( &m_Vacuo)
      );
      serializer.AddDataMember
      (
        "form factor water",
        io::Serialization::GetAgent( &m_Water)
      );

      return serializer;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomGroupTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Vacuo, ISTREAM);
      io::Serialize::Read( m_Mass, ISTREAM);
      io::Serialize::Read( m_SurfaceArea, ISTREAM);
      io::Serialize::Read( m_Water, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &AtomGroupTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Vacuo, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Mass, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SurfaceArea, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Water, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_atom_group_types.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "biol/bcl_biol_atom_group_type_data.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief construct all AtomTypes
    AtomGroupTypes::AtomGroupTypes() :
//                                         Volume   Radius  Bound_Hydrogen,  scattering_length H20, scattering_length D20
//                                         (A^3)    (A)                          (10^-12 cm)             (10^-12 cm)
      H(    AddEnum( "H",   AtomGroupTypeData( GetAtomTypes().H,    5.15, 1.07,  0, -0.3742,  0.6671 ))), //!< Hydrogen
      C(    AddEnum( "C",   AtomGroupTypeData( GetAtomTypes().C,   16.44, 1.58,  0,  0.6651,  0.6651 ))), //!< Carbon
      CH(   AddEnum( "CH",  AtomGroupTypeData( GetAtomTypes().C,   21.59, 1.73,  1,  0.2909,  0.2909 ))), //!< Carbon 1 bound hydrogen
      CH2(  AddEnum( "CH2", AtomGroupTypeData( GetAtomTypes().C,   26.74, 1.85,  2, -0.0833, -0.0833 ))), //!< Carbon 2 bound hydrogens
      CH3(  AddEnum( "CH3", AtomGroupTypeData( GetAtomTypes().C,   31.89, 1.97,  3, -0.4575, -0.4575 ))), //!< Carbon 3 bound hydrogens
      N(    AddEnum( "N",   AtomGroupTypeData( GetAtomTypes().N,    2.49, 0.84,  0,  0.9400,  0.9400 ))), //!< Nitrogen
      NH(   AddEnum( "NH",  AtomGroupTypeData( GetAtomTypes().N,    7.64, 1.22,  1,  0.5658,  1.6071 ))), //!< Nitrogen 1 bound hydrogen
      NH2(  AddEnum( "NH2", AtomGroupTypeData( GetAtomTypes().N,   12.79, 1.45,  2,  0.1916,  2.2742 ))), //!< Nitrogen 2 bound hydrogens
      NH3(  AddEnum( "NH3", AtomGroupTypeData( GetAtomTypes().N,   17.94, 1.62,  3, -0.1826,  2.9413 ))), //!< Nitrogen 3 bound hydrogens
      O(    AddEnum( "O",   AtomGroupTypeData( GetAtomTypes().O,    9.13, 1.30,  0,  0.5804,  0.5804 ))), //!< Oxygen
      OH(   AddEnum( "OH",  AtomGroupTypeData( GetAtomTypes().O,   14.28, 1.50,  1,  0.2062,  1.2475 ))), //!< Oxygen 1 bound hydrogen
      S(    AddEnum( "S",   AtomGroupTypeData( GetAtomTypes().SG,  19.86, 1.68,  0,  0.2847,  0.2847 ))), //!< Sulfur
      SH(   AddEnum( "SH",  AtomGroupTypeData( GetAtomTypes().SG,  25.10, 1.81,  1, -0.0895,  1.2365 ))), //!< Sulfur 1 bound hydrogen
      SE(   AddEnum( "SE",  AtomGroupTypeData( GetAtomTypes().SE,  28.73, 1.90,  0,  0.7971,  0.7971 )))  //!< Selenium
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomGroupTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to the only instance of AtomTypes
    //! @return reference to only instance of AtomTypes
    const AtomGroupTypes &GetAtomGroupTypes()
    {
      return AtomGroupTypes::GetEnums();
    }

    //! @brief map to hold Solvent Volume and Radius values for each heavy atom in each amino acid
    //! @param AA The amino acid types
    //! @param ATOMTYPE The atomtype
    //! @return Group Types
    //! The values are taken from Crysol - a program to evaluate x-ray solution scattering of biological
    //! macromolecules from Atomic Coordinates.  D. Svergun et all 1995
    const AtomGroupType &AtomGroupTypes::GetType( const AAType &AA, const AtomType &ATOMTYPE)
    {
      static storage::Map< AAType, storage::Map< AtomType, AtomGroupType> > map;
      if( map.IsEmpty())
      {
        // Compose Map
        map[ GetAATypes().ALA][ GetAtomTypes().N  ] = GetAtomGroupTypes().N;
        map[ GetAATypes().ALA][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().ALA][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ALA][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().ALA][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().ARG][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().ARG][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().ARG][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ARG][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().ARG][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().ARG][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().ARG][ GetAtomTypes().CD ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().ARG][ GetAtomTypes().NE ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().ARG][ GetAtomTypes().CZ ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ARG][ GetAtomTypes().NH1] = GetAtomGroupTypes().NH2;
        map[ GetAATypes().ARG][ GetAtomTypes().NH2] = GetAtomGroupTypes().NH2;
        map[ GetAATypes().ASN][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().ASN][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().ASN][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ASN][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().ASN][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().ASN][ GetAtomTypes().CG ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ASN][ GetAtomTypes().OD1] = GetAtomGroupTypes().O;
        map[ GetAATypes().ASN][ GetAtomTypes().ND2] = GetAtomGroupTypes().NH2;
        map[ GetAATypes().ASP][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().ASP][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().ASP][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ASP][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().ASP][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().ASP][ GetAtomTypes().CG ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ASP][ GetAtomTypes().OD1] = GetAtomGroupTypes().O;
        map[ GetAATypes().ASP][ GetAtomTypes().OD2] = GetAtomGroupTypes().OH;
        map[ GetAATypes().CYS][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().CYS][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().CYS][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().CYS][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().CYS][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().CYS][ GetAtomTypes().SG ] = GetAtomGroupTypes().SH;
        map[ GetAATypes().GLN][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().GLN][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().GLN][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().GLN][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().GLN][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().GLN][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().GLN][ GetAtomTypes().CD ] = GetAtomGroupTypes().C;
        map[ GetAATypes().GLN][ GetAtomTypes().OE1] = GetAtomGroupTypes().O;
        map[ GetAATypes().GLN][ GetAtomTypes().NE2] = GetAtomGroupTypes().NH2;
        map[ GetAATypes().GLU][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().GLU][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().GLU][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().GLU][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().GLU][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().GLU][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().GLU][ GetAtomTypes().CD ] = GetAtomGroupTypes().C;
        map[ GetAATypes().GLU][ GetAtomTypes().OE1] = GetAtomGroupTypes().O;
        map[ GetAATypes().GLU][ GetAtomTypes().OE2] = GetAtomGroupTypes().OH;
        map[ GetAATypes().GLY][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().GLY][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().GLY][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().GLY][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().HIS][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().HIS][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().HIS][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().HIS][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().HIS][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().HIS][ GetAtomTypes().CG ] = GetAtomGroupTypes().C;
        map[ GetAATypes().HIS][ GetAtomTypes().ND1] = GetAtomGroupTypes().NH;
        map[ GetAATypes().HIS][ GetAtomTypes().CD2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().HIS][ GetAtomTypes().CE1] = GetAtomGroupTypes().CH;
        map[ GetAATypes().HIS][ GetAtomTypes().NE2] = GetAtomGroupTypes().NH;
        map[ GetAATypes().ILE][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().ILE][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().ILE][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ILE][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().ILE][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().ILE][ GetAtomTypes().CG1] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().ILE][ GetAtomTypes().CG2] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().ILE][ GetAtomTypes().CD1] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().LEU][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().LEU][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().LEU][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().LEU][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().LEU][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().LEU][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().LEU][ GetAtomTypes().CD1] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().LEU][ GetAtomTypes().CD2] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().LYS][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().LYS][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().LYS][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().LYS][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().LYS][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().LYS][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().LYS][ GetAtomTypes().CD ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().LYS][ GetAtomTypes().CE ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().LYS][ GetAtomTypes().NZ ] = GetAtomGroupTypes().NH3;
        map[ GetAATypes().MET][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().MET][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().MET][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().MET][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().MET][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().MET][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().MET][ GetAtomTypes().SD ] = GetAtomGroupTypes().S;
        map[ GetAATypes().MET][ GetAtomTypes().CE ] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().MSE][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().MSE][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().MSE][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().MSE][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().MSE][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().MSE][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().MSE][ GetAtomTypes().SE ] = GetAtomGroupTypes().SE;
        map[ GetAATypes().MSE][ GetAtomTypes().CE ] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().PHE][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().PHE][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PHE][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().PHE][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().PHE][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().PHE][ GetAtomTypes().CG ] = GetAtomGroupTypes().C;
        map[ GetAATypes().PHE][ GetAtomTypes().CD1] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PHE][ GetAtomTypes().CD2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PHE][ GetAtomTypes().CE1] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PHE][ GetAtomTypes().CE2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PHE][ GetAtomTypes().CZ ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PRO][ GetAtomTypes().N  ] = GetAtomGroupTypes().N;
        map[ GetAATypes().PRO][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PRO][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().PRO][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().PRO][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().PRO][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().PRO][ GetAtomTypes().CD ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().SER][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().SER][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().SER][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().SER][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().SER][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().SER][ GetAtomTypes().OG ] = GetAtomGroupTypes().OH;
        map[ GetAATypes().THR][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().THR][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().THR][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().THR][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().THR][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().THR][ GetAtomTypes().OG1] = GetAtomGroupTypes().OH;
        map[ GetAATypes().THR][ GetAtomTypes().CG2] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().TRP][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().TRP][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TRP][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().TRP][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().TRP][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().TRP][ GetAtomTypes().CG ] = GetAtomGroupTypes().C;
        map[ GetAATypes().TRP][ GetAtomTypes().CD1] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TRP][ GetAtomTypes().CD2] = GetAtomGroupTypes().C;
        map[ GetAATypes().TRP][ GetAtomTypes().NE1] = GetAtomGroupTypes().NH;
        map[ GetAATypes().TRP][ GetAtomTypes().CE2] = GetAtomGroupTypes().C;
        map[ GetAATypes().TRP][ GetAtomTypes().CE3] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TRP][ GetAtomTypes().CZ2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TRP][ GetAtomTypes().CZ3] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TRP][ GetAtomTypes().CH2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TYR][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().TYR][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TYR][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().TYR][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().TYR][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().TYR][ GetAtomTypes().CG ] = GetAtomGroupTypes().C;
        map[ GetAATypes().TYR][ GetAtomTypes().CD1] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TYR][ GetAtomTypes().CD2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TYR][ GetAtomTypes().CE1] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TYR][ GetAtomTypes().CE2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TYR][ GetAtomTypes().CZ ] = GetAtomGroupTypes().C;
        map[ GetAATypes().TYR][ GetAtomTypes().OH ] = GetAtomGroupTypes().OH;
        map[ GetAATypes().VAL][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().VAL][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().VAL][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().VAL][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().VAL][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().VAL][ GetAtomTypes().CG1] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().VAL][ GetAtomTypes().CG2] = GetAtomGroupTypes().CH3;
      }
      // return group type
      return map.GetValue( AA).GetValue( ATOMTYPE);
    }

    //! @brief map to hold Solvent Volume and Radius values for each heavy atom in each amino acid
    //! @param AA The amino acid types
    //! @param ATOMTYPE The atomtype
    //! @return Group Types
    //! The values are taken from Crysol - a program to evaluate x-ray solution scattering of biological
    //! macromolecules from Atomic Coordinates.  D. Svergun et all 1995
    const std::string &AtomGroupTypes::GetTypeString( const AAType &AA, const AtomType &ATOMTYPE)
    {
      static storage::Map< AAType, storage::Map< AtomType, std::string> > map;
      if( map.IsEmpty())
      {
        // Compose Map
        map[ GetAATypes().ALA][ GetAtomTypes().N  ] = "N";
        map[ GetAATypes().ALA][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().ALA][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().ALA][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().ALA][ GetAtomTypes().CB ] = "CH3";
        map[ GetAATypes().ARG][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().ARG][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().ARG][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().ARG][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().ARG][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().ARG][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().ARG][ GetAtomTypes().CD ] = "CH2";
        map[ GetAATypes().ARG][ GetAtomTypes().NE ] = "NH";
        map[ GetAATypes().ARG][ GetAtomTypes().CZ ] = "C";
        map[ GetAATypes().ARG][ GetAtomTypes().NH1] = "NH2";
        map[ GetAATypes().ARG][ GetAtomTypes().NH2] = "NH2";
        map[ GetAATypes().ASN][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().ASN][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().ASN][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().ASN][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().ASN][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().ASN][ GetAtomTypes().CG ] = "C";
        map[ GetAATypes().ASN][ GetAtomTypes().OD1] = "O";
        map[ GetAATypes().ASN][ GetAtomTypes().ND2] = "NH2";
        map[ GetAATypes().ASP][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().ASP][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().ASP][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().ASP][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().ASP][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().ASP][ GetAtomTypes().CG ] = "C";
        map[ GetAATypes().ASP][ GetAtomTypes().OD1] = "O";
        map[ GetAATypes().ASP][ GetAtomTypes().OD2] = "OH";
        map[ GetAATypes().CYS][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().CYS][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().CYS][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().CYS][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().CYS][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().CYS][ GetAtomTypes().SG ] = "SH";
        map[ GetAATypes().GLN][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().GLN][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().GLN][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().GLN][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().GLN][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().GLN][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().GLN][ GetAtomTypes().CD ] = "C";
        map[ GetAATypes().GLN][ GetAtomTypes().OE1] = "O";
        map[ GetAATypes().GLN][ GetAtomTypes().NE2] = "NH2";
        map[ GetAATypes().GLU][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().GLU][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().GLU][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().GLU][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().GLU][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().GLU][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().GLU][ GetAtomTypes().CD ] = "C";
        map[ GetAATypes().GLU][ GetAtomTypes().OE1] = "O";
        map[ GetAATypes().GLU][ GetAtomTypes().OE2] = "OH";
        map[ GetAATypes().GLY][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().GLY][ GetAtomTypes().CA ] = "CH2";
        map[ GetAATypes().GLY][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().GLY][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().HIS][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().HIS][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().HIS][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().HIS][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().HIS][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().HIS][ GetAtomTypes().CG ] = "C";
        map[ GetAATypes().HIS][ GetAtomTypes().ND1] = "NH";
        map[ GetAATypes().HIS][ GetAtomTypes().CD2] = "CH";
        map[ GetAATypes().HIS][ GetAtomTypes().CE1] = "CH";
        map[ GetAATypes().HIS][ GetAtomTypes().NE2] = "NH";
        map[ GetAATypes().ILE][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().ILE][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().ILE][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().ILE][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().ILE][ GetAtomTypes().CB ] = "CH";
        map[ GetAATypes().ILE][ GetAtomTypes().CG1] = "CH2";
        map[ GetAATypes().ILE][ GetAtomTypes().CG2] = "CH3";
        map[ GetAATypes().ILE][ GetAtomTypes().CD1] = "CH3";
        map[ GetAATypes().LEU][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().LEU][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().LEU][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().LEU][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().LEU][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().LEU][ GetAtomTypes().CG ] = "CH";
        map[ GetAATypes().LEU][ GetAtomTypes().CD1] = "CH3";
        map[ GetAATypes().LEU][ GetAtomTypes().CD2] = "CH3";
        map[ GetAATypes().LYS][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().LYS][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().LYS][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().LYS][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().LYS][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().LYS][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().LYS][ GetAtomTypes().CD ] = "CH2";
        map[ GetAATypes().LYS][ GetAtomTypes().CE ] = "CH2";
        map[ GetAATypes().LYS][ GetAtomTypes().NZ ] = "NH3";
        map[ GetAATypes().MET][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().MET][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().MET][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().MET][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().MET][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().MET][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().MET][ GetAtomTypes().SD ] = "S";
        map[ GetAATypes().MET][ GetAtomTypes().CE ] = "CH3";
        map[ GetAATypes().MSE][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().MSE][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().MSE][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().MSE][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().MSE][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().MSE][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().MSE][ GetAtomTypes().SE ] = "SE";
        map[ GetAATypes().MSE][ GetAtomTypes().CE ] = "CH3";
        map[ GetAATypes().PHE][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().PHE][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().PHE][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().PHE][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().PHE][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().PHE][ GetAtomTypes().CG ] = "C";
        map[ GetAATypes().PHE][ GetAtomTypes().CD1] = "CH";
        map[ GetAATypes().PHE][ GetAtomTypes().CD2] = "CH";
        map[ GetAATypes().PHE][ GetAtomTypes().CE1] = "CH";
        map[ GetAATypes().PHE][ GetAtomTypes().CE2] = "CH";
        map[ GetAATypes().PHE][ GetAtomTypes().CZ ] = "CH";
        map[ GetAATypes().PRO][ GetAtomTypes().N  ] = "N";
        map[ GetAATypes().PRO][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().PRO][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().PRO][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().PRO][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().PRO][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().PRO][ GetAtomTypes().CD ] = "CH2";
        map[ GetAATypes().SER][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().SER][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().SER][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().SER][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().SER][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().SER][ GetAtomTypes().OG ] = "OH";
        map[ GetAATypes().THR][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().THR][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().THR][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().THR][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().THR][ GetAtomTypes().CB ] = "CH";
        map[ GetAATypes().THR][ GetAtomTypes().OG1] = "OH";
        map[ GetAATypes().THR][ GetAtomTypes().CG2] = "CH3";
        map[ GetAATypes().TRP][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().TRP][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().TRP][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().TRP][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().TRP][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().TRP][ GetAtomTypes().CG ] = "C";
        map[ GetAATypes().TRP][ GetAtomTypes().CD1] = "CH";
        map[ GetAATypes().TRP][ GetAtomTypes().CD2] = "C";
        map[ GetAATypes().TRP][ GetAtomTypes().NE1] = "NH";
        map[ GetAATypes().TRP][ GetAtomTypes().CE2] = "C";
        map[ GetAATypes().TRP][ GetAtomTypes().CE3] = "CH";
        map[ GetAATypes().TRP][ GetAtomTypes().CZ2] = "CH";
        map[ GetAATypes().TRP][ GetAtomTypes().CZ3] = "CH";
        map[ GetAATypes().TRP][ GetAtomTypes().CH2] = "CH";
        map[ GetAATypes().TYR][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().TYR][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().TYR][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().TYR][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().TYR][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().TYR][ GetAtomTypes().CG ] = "C";
        map[ GetAATypes().TYR][ GetAtomTypes().CD1] = "CH";
        map[ GetAATypes().TYR][ GetAtomTypes().CD2] = "CH";
        map[ GetAATypes().TYR][ GetAtomTypes().CE1] = "CH";
        map[ GetAATypes().TYR][ GetAtomTypes().CE2] = "CH";
        map[ GetAATypes().TYR][ GetAtomTypes().CZ ] = "C";
        map[ GetAATypes().TYR][ GetAtomTypes().OH ] = "OH";
        map[ GetAATypes().VAL][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().VAL][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().VAL][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().VAL][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().VAL][ GetAtomTypes().CB ] = "CH";
        map[ GetAATypes().VAL][ GetAtomTypes().CG1] = "CH3";
        map[ GetAATypes().VAL][ GetAtomTypes().CG2] = "CH3";
      }
      // return group type
      return map.GetValue( AA).GetValue( ATOMTYPE);
    }

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< biol::AtomGroupTypeData, biol::AtomGroupTypes>;

  } // namespace util
} // namespace bcl
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
#include "biol/bcl_biol_atom_type_data.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AtomTypeData::s_Instance( GetObjectInstances().AddInstance( new AtomTypeData()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct undefined atom type
    AtomTypeData::AtomTypeData() :
      m_AtomName( ""),
      m_ElementType( util::GetUndefined< chemistry::ElementType>()),
      m_InBackBone( false),
      m_InSideChain( false),
      m_HelixCoordinates( util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>()),
      m_StrandCoordinates( util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>()),
      m_BondLengths(),
      m_Connections()
    {
    }

    //! @brief construct undefined atom type from given information
    //! @param ATOM_NAME name of the atom
    //! @param ELEMENT_TYPE element type of the atom
    //! @param IN_BACKBONE if atom type is in backbone
    //! @param IN_SIDE_CHAIN true if atom type is in side chain
    //! @param HELIX_CYLINDER_COORDINATES position on the cylinder for helix
    //! @param STRAND_CYLINDER_COORDINATES position on the cylinder for strand
    AtomTypeData::AtomTypeData
    (
      const std::string &ATOM_NAME,
      const chemistry::ElementType &ELEMENT_TYPE,
      const bool IN_BACKBONE,
      const bool IN_SIDE_CHAIN,
      const coord::CylinderCoordinates &HELIX_CYLINDER_COORDINATES,
      const coord::CylinderCoordinates &STRAND_CYLINDER_COORDINATES
    ) :
      m_AtomName( ATOM_NAME),
      m_ElementType( ELEMENT_TYPE),
      m_InBackBone( IN_BACKBONE),
      m_InSideChain( IN_SIDE_CHAIN),
      m_HelixCoordinates( HELIX_CYLINDER_COORDINATES),
      m_StrandCoordinates( STRAND_CYLINDER_COORDINATES),
      m_BondLengths(),
      m_Connections()
    {
    }

    //! @brief virtual copy constructor
    AtomTypeData *AtomTypeData::Clone() const
    {
      return new AtomTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return Cylinder Coordinates for the specified SS_TYPE
    //! @param SS_TYPE SSType of interest
    //! @return Cylinder Coordinates for the specified SS_TYPE
    const coord::CylinderCoordinates &AtomTypeData::GetCoordinates( const SSType &SS_TYPE) const
    {
      static coord::CylinderCoordinates s_coordinates;

      if( SS_TYPE == GetSSTypes().HELIX)
      {
        return m_HelixCoordinates;
      }
      if( SS_TYPE == GetSSTypes().STRAND)
      {
        return m_StrandCoordinates;
      }
      BCL_Exit( "No cylindrical coordinates are stored for the provided sstype " + util::Format()( SS_TYPE), -1);

      return s_coordinates;
    }

    //! @brief return Cylinder Coordinates for the specified SS_TYPE
    //! @param SS_TYPE SSType of interest
    //! @return Cylinder Coordinates for the specified SS_TYPE
    const storage::Set< AtomType> &AtomTypeData::GetConnections() const
    {
      return m_Connections;
    }

    //! @brief return atom types that this atom type will connect to with a double bond, if present in the AA
    //! @return atom types that this atom type will connect to with a double bond, if present in the AA
    const storage::Set< AtomType> &AtomTypeData::GetDoubleBondConnections() const
    {
      return m_DoubleBondConnections;
    }

    //! @brief gets the bond length between this atom type and the passed atom type
    //! @param ATOM_TYPE atom type bound to
    //! @return the bond length between this atom type and the passed atom type
    double AtomTypeData::GetBondLength( const AtomType &ATOM_TYPE) const
    {
      // search the map for the atom type
      const storage::Map< AtomType, double>::const_iterator find_itr( m_BondLengths.Find( ATOM_TYPE));

      // return an undefined double if it wasn't found, otherwise return the length
      return find_itr == m_BondLengths.End() ? util::GetUndefined< double>() : find_itr->second;
    }

    //! @brief sets the bond length to another atom type
    //! @param ATOM_TYPE atom type bound to
    //! @param LENGTH average bond length
    void AtomTypeData::SetBondLength( const AtomType &ATOM_TYPE, const double LENGTH) const
    {
      m_BondLengths[ ATOM_TYPE] = LENGTH;
    }

    //! @brief sets the connections that the atom always makes, if available
    //! @param CONNECTIONS set of atom types that this atom type always connects to
    void AtomTypeData::SetConnections( const storage::Set< AtomType> &CONNECTIONS) const
    {
      m_Connections = CONNECTIONS;
    }

    //! @brief sets the double bond connections that the atom always makes, if available
    //! @param CONNECTIONS set of atom types that this atom type always makes double bonds to
    void AtomTypeData::SetDoubleBondConnections( const storage::Set< AtomType> &CONNECTIONS) const
    {
      m_DoubleBondConnections = CONNECTIONS;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_AtomName, ISTREAM);
      io::Serialize::Read( m_ElementType, ISTREAM);
      io::Serialize::Read( m_InBackBone, ISTREAM);
      io::Serialize::Read( m_HelixCoordinates, ISTREAM);
      io::Serialize::Read( m_StrandCoordinates, ISTREAM);
      io::Serialize::Read( m_BondLengths, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &AtomTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_AtomName, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ElementType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_InBackBone, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HelixCoordinates, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_StrandCoordinates, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BondLengths, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_atom_types.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief construct all AtomTypes
    AtomTypes::AtomTypes() :
//                                                 ElementType                              in BB  in SC, HELIX                                                                   STRAND
      N(    AddEnum( "N",    AtomTypeData( " N  ", chemistry::GetElementTypes().e_Nitrogen, true,  false, coord::CylinderCoordinates( -0.245, 1.520, -27.600 / 180 * math::g_Pi), coord::CylinderCoordinates( 0.477, 0.391,  0.998)))), //!< Nitrogen from the peptide bond
      CA(   AddEnum( "CA",   AtomTypeData( " CA ", chemistry::GetElementTypes().e_Carbon,   true,  false, coord::CylinderCoordinates(  0.661, 2.260,   0.000 / 180 * math::g_Pi), coord::CylinderCoordinates( 1.746, 0.778, -0.377)))), //!< Carbon alpha backbone
      C(    AddEnum( "C",    AtomTypeData( " C  ", chemistry::GetElementTypes().e_Carbon,   true,  false, coord::CylinderCoordinates(  1.747, 1.678,  26.335 / 180 * math::g_Pi), coord::CylinderCoordinates( 2.949, 0.534,  1.126)))), //!< Carbon from the carboxyl group
      O(    AddEnum( "O",    AtomTypeData( " O  ", chemistry::GetElementTypes().e_Oxygen,   true,  false, coord::CylinderCoordinates(  2.922, 1.950,  20.250 / 180 * math::g_Pi), coord::CylinderCoordinates( 3.001, 1.720,  1.424)))), //!< Oxygen from the carboxyl group
      CB(   AddEnum( "CB",   AtomTypeData( " CB ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates( -0.151, 3.240,  18.680 / 180 * math::g_Pi), coord::CylinderCoordinates( 1.756, 2.275, -0.106)))), //!< Carbon beta first side chain atom
      CG(   AddEnum( "CG",   AtomTypeData( " CG ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon gamma - second side chain atom
      CG1(  AddEnum( "CG1",  AtomTypeData( " CG1", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon gamma - second side chain atom 1 for two second side chain atoms
      CG2(  AddEnum( "CG2",  AtomTypeData( " CG2", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon gamma - second side chain atom 2 for two second side chain atoms
      CD(   AddEnum( "CD",   AtomTypeData( " CD ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon delta - third side chain atom
      CD1(  AddEnum( "CD1",  AtomTypeData( " CD1", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon delta - third side chain atom 1 for two third side chain atoms
      CD2(  AddEnum( "CD2",  AtomTypeData( " CD2", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon delta - third side chain atom 2 for two third side chain atoms
      CE(   AddEnum( "CE",   AtomTypeData( " CE ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon epsilon - fourth side chain atom
      CE1(  AddEnum( "CE1",  AtomTypeData( " CE1", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon epsilon - fourth side chain atom 1 for two or three fourth side chain atoms
      CE2(  AddEnum( "CE2",  AtomTypeData( " CE2", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon epsilon - fourth side chain atom 2 for two or three fourth side chain atoms
      CE3(  AddEnum( "CE3",  AtomTypeData( " CE3", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon epsilon - fourth side chain atom 3 for two or three fourth side chain atoms
      CZ(   AddEnum( "CZ",   AtomTypeData( " CZ ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon zeta - fifth side chain atom
      CZ2(  AddEnum( "CZ2",  AtomTypeData( " CZ2", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon zeta - fifth side chain atom 1 for two fifth side chain atoms
      CZ3(  AddEnum( "CZ3",  AtomTypeData( " CZ3", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon zeta - fifth side chain atom 2 for two fifth side chain atoms
      CH2(  AddEnum( "CH2",  AtomTypeData( " CH2", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon eta - sixth side chain aom as in TRP
      ND1(  AddEnum( "ND1",  AtomTypeData( " ND1", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen delta - third position
      ND2(  AddEnum( "ND2",  AtomTypeData( " ND2", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen delta - third position
      NE(   AddEnum( "NE",   AtomTypeData( " NE ", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen epsilon - fourth position nitrogen as in ARG
      NE1(  AddEnum( "NE1",  AtomTypeData( " NE1", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen epsilon - fourth position nitrogen as in TRP
      NE2(  AddEnum( "NE2",  AtomTypeData( " NE2", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen epsilon - fourth position nitrogen as in GLN or HIS
      NZ(   AddEnum( "NZ",   AtomTypeData( " NZ ", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen zeta - fifth side chain atom as in LYS
      NH1(  AddEnum( "NH1",  AtomTypeData( " NH1", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen eta 1 - sixth side chain atom as in ARG
      NH2(  AddEnum( "NH2",  AtomTypeData( " NH2", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen eta 2 - sixth side chain atom as in ARG
      OD1(  AddEnum( "OD1",  AtomTypeData( " OD1", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen delta - third side chain atom as in ASP or ASN
      OD2(  AddEnum( "OD2",  AtomTypeData( " OD2", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen delta - third side chain atom as in ASP
      OG(   AddEnum( "OG",   AtomTypeData( " OG ", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen gamma - second side chain atom as in SER
      OG1(  AddEnum( "OG1",  AtomTypeData( " OG1", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen gamma - second side chain atom as in THR
      OE1(  AddEnum( "OE1",  AtomTypeData( " OE1", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen epsilon - fourth side chain atom as in GLN
      OE2(  AddEnum( "OE2",  AtomTypeData( " OE2", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen epsilon - fourth side chain atom as in GLN or GLU
      OH(   AddEnum( "OH",   AtomTypeData( " OH ", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen eta - sixth side chain atom as in TYR
      SD(   AddEnum( "SD",   AtomTypeData( " SD ", chemistry::GetElementTypes().e_Sulfur,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Suflur on delta carbon - like in methionin
      SE(   AddEnum( "SE",   AtomTypeData( " SE ", chemistry::GetElementTypes().e_Selenium, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Selenium on delta carbon - like in selenomethionine
      SG(   AddEnum( "SG",   AtomTypeData( " SG ", chemistry::GetElementTypes().e_Sulfur,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Sulfur on gamma carbon - like in cystein
      H(    AddEnum( "H",    AtomTypeData( " H  ", chemistry::GetElementTypes().e_Hydrogen, false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on the backbone nitrogen atom
      HA(   AddEnum( "HA",   AtomTypeData( " HA ", chemistry::GetElementTypes().e_Hydrogen, false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CA alpha Carbon
      HA2(  AddEnum( "HA2",  AtomTypeData( " HA2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates( -0.151, 3.240,  18.680 / 180 * math::g_Pi), coord::CylinderCoordinates( 1.756, 2.275, -0.106)))), //!< Hydrogen on CA for GLY in CB position
      HA3(  AddEnum( "HA3",  AtomTypeData( " HA3", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CA for GLY in HA position
      HB(   AddEnum( "HB",   AtomTypeData( " HB ", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CB for THR, ILE or VAL
      HB1(  AddEnum( "HB1",  AtomTypeData( " HB1", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CB for ALA
      HB2(  AddEnum( "HB2",  AtomTypeData( " HB2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CB
      HB3(  AddEnum( "HB3",  AtomTypeData( " HB3", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CB
      HG(   AddEnum( "HG",   AtomTypeData( " HG ", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen an CG as in LEU, CYS, SER
      HG1(  AddEnum( "HG1",  AtomTypeData( " HG1", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 an CG as in THR
      HG2(  AddEnum( "HG2",  AtomTypeData( " HG2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CG
      HG3(  AddEnum( "HG3",  AtomTypeData( " HG3", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CG
      HG11( AddEnum( "HG11", AtomTypeData( "HG11", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CG1 as in VAL
      HG12( AddEnum( "HG12", AtomTypeData( "HG12", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CG2 as in ILE, VAL
      HG13( AddEnum( "HG13", AtomTypeData( "HG13", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CG1 as in ILE, VAL
      HG21( AddEnum( "HG21", AtomTypeData( "HG21", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CG2 as in ILE, VAL, THR
      HG22( AddEnum( "HG22", AtomTypeData( "HG22", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CG2 as in ILE, VAL, THR
      HG23( AddEnum( "HG23", AtomTypeData( "HG23", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CG2 as in ILE, VAL, THR
      HD1(  AddEnum( "HD1",  AtomTypeData( " HD1", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CD
      HD2(  AddEnum( "HD2",  AtomTypeData( " HD2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CD
      HD3(  AddEnum( "HD3",  AtomTypeData( " HD3", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CD
      HD11( AddEnum( "HD11", AtomTypeData( "HD11", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CD1
      HD12( AddEnum( "HD12", AtomTypeData( "HD12", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CD1
      HD13( AddEnum( "HD13", AtomTypeData( "HD13", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CD1
      HD21( AddEnum( "HD21", AtomTypeData( "HD21", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CD2
      HD22( AddEnum( "HD22", AtomTypeData( "HD22", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CD2
      HD23( AddEnum( "HD23", AtomTypeData( "HD23", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CD2
      HE(   AddEnum( "HE",   AtomTypeData( " HE ", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CE as in ARG
      HE1(  AddEnum( "HE1",  AtomTypeData( " HE1", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CE
      HE2(  AddEnum( "HE2",  AtomTypeData( " HE2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CE
      HE3(  AddEnum( "HE3",  AtomTypeData( " HE3", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CE
      HE21( AddEnum( "HE21", AtomTypeData( "HE21", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CE2 as in GLN
      HE22( AddEnum( "HE22", AtomTypeData( "HE22", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CE2 as in GLN
      HZ(   AddEnum( "HZ",   AtomTypeData( " HZ ", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CZ as in PHE
      HZ1(  AddEnum( "HZ1",  AtomTypeData( " HZ1", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on Z1 as in LYS (on Nitrogen NZ)
      HZ2(  AddEnum( "HZ2",  AtomTypeData( " HZ2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on Z as in LYS (on NZ), TRP (on CZ2)
      HZ3(  AddEnum( "HZ3",  AtomTypeData( " HZ3", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on Z as in LYS (on NZ), TRP (on CZ3)
      HH(   AddEnum( "HH",   AtomTypeData( " HH ", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CH as in TYR
      HH2(  AddEnum( "HH2",  AtomTypeData( " HH2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CH2 as in TRP
      HH11( AddEnum( "HH11", AtomTypeData( "HH11", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on NH1 as in ARG
      HH12( AddEnum( "HH12", AtomTypeData( "HH12", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on NH1 as in ARG
      HH21( AddEnum( "HH21", AtomTypeData( "HH21", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on NH2 as in ARG
      HH22( AddEnum( "HH22", AtomTypeData( "HH22", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on NH2 as in ARG
      // terminal amine
      H1(   AddEnum( "H1",   AtomTypeData( " H1 ", chemistry::GetElementTypes().e_Hydrogen, false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< leaving hydrogen, connected to backbone Nitrogen
      H2(   AddEnum( "H2",   AtomTypeData( " H2 ", chemistry::GetElementTypes().e_Hydrogen, false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< leaving hydrogen, connected to backbone Nitrogen
      H3(   AddEnum( "H3",   AtomTypeData( " H3 ", chemistry::GetElementTypes().e_Hydrogen, false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< leaving hydrogen, connected to backbone Nitrogen
      // terminal carboxylic acid
      HXT(  AddEnum( "HXT",  AtomTypeData( " HXT", chemistry::GetElementTypes().e_Hydrogen, false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen leaving hydrogen, connected to C if there is no peptide bond
      OXT(  AddEnum( "OXT",  AtomTypeData( " OXT", chemistry::GetElementTypes().e_Oxygen,   false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< leaving oxygen, connected to C if there is no peptide bond

      C2(   AddEnum( "C2",   AtomTypeData( " C2 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mttsl C next to N-O nearest to linker
      C3(   AddEnum( "C3",   AtomTypeData( " C3 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C in "pentane" ring connected to linker
      C4(   AddEnum( "C4",   AtomTypeData( " C4 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C in pentane ring connected to H and connected to C3 by double bond
      C5(   AddEnum( "C5",   AtomTypeData( " C5 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C next to N-O nearest to linker on side with C=C-H and opposite C2
      C6(   AddEnum( "C6",   AtomTypeData( " C6 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C above C5 when looking at ring from linker to N-O and oriented C=C-H
      C7(   AddEnum( "C7",   AtomTypeData( " C7 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C below C5 when looking at ring from linker to N-O and oriented C=C-H
      C8(   AddEnum( "C8",   AtomTypeData( " C8 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C above C2 when looking at ring from linker to N-O and oriented C=C-H
      C9(   AddEnum( "C9",   AtomTypeData( " C9 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C below C2 when looking at ring from linker to N-O and oriented C=C-H
      N1(   AddEnum( "N1",   AtomTypeData( " N1 ", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl N in ring creating N-O nitroxide moiety
      O1(   AddEnum( "O1",   AtomTypeData( " O1 ", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     ))))  //!< mtssl O in ring creating N-O nitroxide moiety

    {
      // set bond lengths for backbone atom types (taken from rosetta param files)
        N->SetBondLength(   C, 1.329);
        N->SetBondLength(  CA, 1.458);
        N->SetBondLength(   H, 1.03);
       CA->SetBondLength(  HA, 1.09);
       CA->SetBondLength( HA2, 1.09);
       CA->SetBondLength( HA3, 1.09);
       CA->SetBondLength(   N, 1.458);
       CA->SetBondLength(   C, 1.523);
       CA->SetBondLength(  CB, 1.533);
        C->SetBondLength(   N, 1.329);
        C->SetBondLength(   O, 1.231);
        C->SetBondLength(  CA, 1.523);
        O->SetBondLength(   C, 1.231);
       CB->SetBondLength(  CA, 1.533);
        H->SetBondLength(   N, 1.03);
       HA->SetBondLength(  CA, 1.09);
      HA2->SetBondLength(  CA, 1.09);
      HA3->SetBondLength(  CA, 1.09);

      // set obligatory connections for all simple atom types (e.g. connections that will always be made provided that
      // both Atom types are present)
      // The only non-mandatory connection is N <> CD, which only occurs in proline
      // This list was generated (mostly) automatically from the rosetta params files, after verification that all
      // connections other than the proline N <> CD obligatory
      C->SetConnections( storage::Set< AtomType>::Create( CA, O, OXT));
      C->SetDoubleBondConnections( storage::Set< AtomType>::Create( O));
      CA->SetConnections( storage::Set< AtomType>::Create( C, CB, HA, HA2, HA3, N));
      CB->SetConnections( storage::Set< AtomType>::Create( CA, CG, CG1, CG2, HB, HB1, HB2, HB3, OG, OG1, SG));
      CD->SetConnections( storage::Set< AtomType>::Create( CE, CG, HD1, HD2, HD3, NE, NE2, OE1, OE2));
      CD->SetDoubleBondConnections( storage::Set< AtomType>::Create( OE1));
      CD1->SetConnections( storage::Set< AtomType>::Create( CE1, CG, CG1, HD1, HD11, HD12, HD13, NE1));
      CD1->SetDoubleBondConnections( storage::Set< AtomType>::Create( CG));
      CD2->SetConnections( storage::Set< AtomType>::Create( CE2, CE3, CG, HD2, HD21, HD22, HD23, NE2));
      CD2->SetDoubleBondConnections( storage::Set< AtomType>::Create( CE2)); // HIS, TRP, PHE, TYR
      CE->SetConnections( storage::Set< AtomType>::Create( CD, HE1, HE2, HE3, NZ, SD, SE, C3));
      CE1->SetConnections( storage::Set< AtomType>::Create( CD1, CZ, HE1, ND1, NE2));
      CE1->SetDoubleBondConnections( storage::Set< AtomType>::Create( ND1, CD1, CZ));
      CE2->SetConnections( storage::Set< AtomType>::Create( CD2, CZ, CZ2, HE2, NE1));
      CE2->SetDoubleBondConnections( storage::Set< AtomType>::Create( CZ, CD2)); // PHE, TYR
      CE3->SetConnections( storage::Set< AtomType>::Create( CD2, CZ3, HE3));
      CE3->SetDoubleBondConnections( storage::Set< AtomType>::Create( CZ3));
      CG->SetConnections( storage::Set< AtomType>::Create( CB, CD, CD1, CD2, HG, HG1, HG2, HG3, ND1, ND2, OD1, OD2, SD, SE));
      CG->SetDoubleBondConnections( storage::Set< AtomType>::Create( OD1, CD1));
      CG1->SetConnections( storage::Set< AtomType>::Create( CB, CD1, HG11, HG12, HG13));
      CG2->SetConnections( storage::Set< AtomType>::Create( CB, HG21, HG22, HG23));
      CH2->SetConnections( storage::Set< AtomType>::Create( CZ2, CZ3, HH2));
      CH2->SetDoubleBondConnections( storage::Set< AtomType>::Create( CZ2));
      CZ->SetConnections( storage::Set< AtomType>::Create( CE1, CE2, HZ, NE, NH1, NH2, OH));
      CZ->SetDoubleBondConnections( storage::Set< AtomType>::Create( NH1, CE1)); // ARG / PHE, TYR
      CZ2->SetConnections( storage::Set< AtomType>::Create( CE2, CH2, HZ2));
      CZ2->SetDoubleBondConnections( storage::Set< AtomType>::Create( CH2));
      CZ3->SetConnections( storage::Set< AtomType>::Create( CE3, CH2, HZ3));
      CZ3->SetDoubleBondConnections( storage::Set< AtomType>::Create( CE3));
      H->SetConnections( storage::Set< AtomType>::Create( N));
      HA->SetConnections( storage::Set< AtomType>::Create( CA));
      HA2->SetConnections( storage::Set< AtomType>::Create( CA));
      HA3->SetConnections( storage::Set< AtomType>::Create( CA));
      HB->SetConnections( storage::Set< AtomType>::Create( CB));
      HB1->SetConnections( storage::Set< AtomType>::Create( CB));
      HB2->SetConnections( storage::Set< AtomType>::Create( CB));
      HB3->SetConnections( storage::Set< AtomType>::Create( CB));

      // Additional double bonds that are AA-dependent:
      // CG = CD2, not CD1 (Histidine only)

      // Connection ND1 <> HD1 is only present on HIS at pH < 6, but evidently it was decided that pH < 6 is the most
      // important since HIS was given an HD1. The reasons behind this decision should be investigated and documented
      HD1->SetConnections( storage::Set< AtomType>::Create( CD, CD1, ND1));
      HD11->SetConnections( storage::Set< AtomType>::Create( CD1));
      HD12->SetConnections( storage::Set< AtomType>::Create( CD1));
      HD13->SetConnections( storage::Set< AtomType>::Create( CD1));
      HD2->SetConnections( storage::Set< AtomType>::Create( CD, CD2, OD2));
      HD21->SetConnections( storage::Set< AtomType>::Create( CD2, ND2));
      HD22->SetConnections( storage::Set< AtomType>::Create( CD2, ND2));
      HD23->SetConnections( storage::Set< AtomType>::Create( CD2));
      HD3->SetConnections( storage::Set< AtomType>::Create( CD));
      HE->SetConnections( storage::Set< AtomType>::Create( NE));
      HE1->SetConnections( storage::Set< AtomType>::Create( CE, CE1, NE1));
      HE2->SetConnections( storage::Set< AtomType>::Create( CE, CE2, NE2, OE2));
      HE21->SetConnections( storage::Set< AtomType>::Create( NE2));
      HE22->SetConnections( storage::Set< AtomType>::Create( NE2));
      HE3->SetConnections( storage::Set< AtomType>::Create( CE, CE3));
      HG->SetConnections( storage::Set< AtomType>::Create( CG, OG, SG));
      HG1->SetConnections( storage::Set< AtomType>::Create( CG, OG1));
      HG11->SetConnections( storage::Set< AtomType>::Create( CG1));
      HG12->SetConnections( storage::Set< AtomType>::Create( CG1));
      HG13->SetConnections( storage::Set< AtomType>::Create( CG1));
      HG2->SetConnections( storage::Set< AtomType>::Create( CG));
      HG21->SetConnections( storage::Set< AtomType>::Create( CG2));
      HG22->SetConnections( storage::Set< AtomType>::Create( CG2));
      HG23->SetConnections( storage::Set< AtomType>::Create( CG2));
      HG3->SetConnections( storage::Set< AtomType>::Create( CG));
      HH->SetConnections( storage::Set< AtomType>::Create( OH));
      HH11->SetConnections( storage::Set< AtomType>::Create( NH1));
      HH12->SetConnections( storage::Set< AtomType>::Create( NH1));
      HH2->SetConnections( storage::Set< AtomType>::Create( CH2));
      HH21->SetConnections( storage::Set< AtomType>::Create( NH2));
      HH22->SetConnections( storage::Set< AtomType>::Create( NH2));
      HZ->SetConnections( storage::Set< AtomType>::Create( CZ));
      HZ1->SetConnections( storage::Set< AtomType>::Create( NZ));
      HZ2->SetConnections( storage::Set< AtomType>::Create( CZ2, NZ));
      HZ3->SetConnections( storage::Set< AtomType>::Create( CZ3, NZ));
      N->SetConnections( storage::Set< AtomType>::Create( CA, H, H1, H2, H3));
      ND1->SetConnections( storage::Set< AtomType>::Create( CE1, CG, HD1));
      ND1->SetDoubleBondConnections( storage::Set< AtomType>::Create( CE1));
      ND2->SetConnections( storage::Set< AtomType>::Create( CG, HD21, HD22));
      NE->SetConnections( storage::Set< AtomType>::Create( CD, CZ, HE));
      NE1->SetConnections( storage::Set< AtomType>::Create( CD1, CE2, HE1));
      NE2->SetConnections( storage::Set< AtomType>::Create( CD, CD2, CE1, HE2, HE21, HE22));
      NH1->SetConnections( storage::Set< AtomType>::Create( CZ, HH11, HH12));
      NH1->SetDoubleBondConnections( storage::Set< AtomType>::Create( CZ)); // ARG
      NH2->SetConnections( storage::Set< AtomType>::Create( CZ, HH21, HH22));
      NZ->SetConnections( storage::Set< AtomType>::Create( CE, HZ1, HZ2, HZ3));
      O->SetConnections( storage::Set< AtomType>::Create( C));
      O->SetDoubleBondConnections( storage::Set< AtomType>::Create( C)); // Backbone O
      OD1->SetConnections( storage::Set< AtomType>::Create( CG));
      OD1->SetDoubleBondConnections( storage::Set< AtomType>::Create( CG)); // ASP, ASN
      OD2->SetConnections( storage::Set< AtomType>::Create( CG, HD2));
      OE1->SetConnections( storage::Set< AtomType>::Create( CD));
      OE1->SetDoubleBondConnections( storage::Set< AtomType>::Create( CD)); // GLN, GLU
      OE2->SetConnections( storage::Set< AtomType>::Create( CD, HE2));
      OG->SetConnections( storage::Set< AtomType>::Create( CB, HG));
      OG1->SetConnections( storage::Set< AtomType>::Create( CB, HG1));
      OH->SetConnections( storage::Set< AtomType>::Create( CZ, HH));
      SD->SetConnections( storage::Set< AtomType>::Create( CE, CG, SG));
      SE->SetConnections( storage::Set< AtomType>::Create( CG, CE));
      SG->SetConnections( storage::Set< AtomType>::Create( CB, HG, SD));

      // terminal amine
      H1->SetConnections( storage::Set< AtomType>::Create( N));
      H2->SetConnections( storage::Set< AtomType>::Create( N));
      H3->SetConnections( storage::Set< AtomType>::Create( N));

      // terminal carboxylic acid
      HXT->SetConnections( storage::Set< AtomType>::Create( OXT));
      OXT->SetConnections( storage::Set< AtomType>::Create( C, HXT));

      // methanesulfonothioate - specific
      C2->SetConnections( storage::Set< AtomType>::Create( C3, C8, C9, N1));
      C3->SetConnections( storage::Set< AtomType>::Create( C2, C4, CE));
      C3->SetDoubleBondConnections( storage::Set< AtomType>::Create( C4));
      C4->SetConnections( storage::Set< AtomType>::Create( C3, C5));
      C4->SetDoubleBondConnections( storage::Set< AtomType>::Create( C3));
      C5->SetConnections( storage::Set< AtomType>::Create( C4, C6, C7, N1));
      C6->SetConnections( storage::Set< AtomType>::Create( C5));
      C7->SetConnections( storage::Set< AtomType>::Create( C5));
      C8->SetConnections( storage::Set< AtomType>::Create( C2));
      C9->SetConnections( storage::Set< AtomType>::Create( C2));
      N1->SetConnections( storage::Set< AtomType>::Create( C2, C5, O1));
    }

    //! @brief conversion to a string from a Subset
    //! @param SUBSET the subset to get a string for
    //! @return a string representing that subset
    const std::string &AtomTypes::GetSubsetName( const Subset &SUBSET)
    {
      static const std::string s_descriptors[] =
      {
        "All",
        "Backbone",
        "Side-chain"
      };
      return s_descriptors[ size_t( SUBSET)];
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return command line flag for defining the atoms used in the calculation
    //! @return command line flag for defining the the atoms used in the calculation
    util::ShPtr< command::FlagInterface> &AtomTypes::GetFlagAtomTypes()
    {
      // initialize atom types
      storage::Vector< std::string> atom_types
      (
        storage::Vector< std::string>::Create
        (
          GetSubsetName( e_All),
          GetSubsetName( e_Backbone),
          GetSubsetName( e_SideChain)
        )
      );
      atom_types.Append( GetAtomTypes().GetEnumStrings());

      // initialize flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "atoms", "list of atoms to be included in quality calculation, CA by default",
           command::Parameter
          (
            "atom",
            "any atom from the list (backbone and side chain atoms)",
            command::ParameterCheckAllowed( atom_types),
            GetAtomTypes().CA.GetName()
          ),
          0,
          GetAtomTypes().GetEnumCount()
        )
      );

      // end
      return s_flag;
    }

    //! @brief get the set of atoms defined by the flag
    //! @return the set of atoms defined by the flag
    storage::Set< AtomType> AtomTypes::GetCommandLineAtoms()
    {
      // check for "All"
      const storage::Set< std::string> params( GetFlagAtomTypes()->GetObjectSet< std::string>());
      if( params.Contains( GetSubsetName( e_All)))
      {
        // return all atom types
        return storage::Set< AtomType>( GetAtomTypes().Begin(), GetAtomTypes().End());
      }
      else if( params.Contains( GetSubsetName( e_Backbone)))
      {
        // return backbone heavy atom types
        return GetAtomTypes().GetBackBoneAtomTypes();
      }
      else if( params.Contains( GetSubsetName( e_SideChain)))
      {
        // return side chain heavy atom types
        return GetAtomTypes().GetSideChainAtomTypes();
      }

      // get specified atom types
      storage::Set< AtomType> atom_types( GetFlagAtomTypes()->GetObjectSet< AtomType>());

      // add CA if empty
      if( atom_types.IsEmpty())
      {
        atom_types.Insert( GetAtomTypes().CA);
      }

      // end
      return atom_types;
    }

    //! @brief return StorageVector of Backbone Atom Types
    //! @return StorageVector of Backbone Atom Types
    const storage::Set< AtomType> &AtomTypes::GetBackBoneAtomTypes() const
    {
      // initialize static set of AtomTypes from first enum until CB
      static const storage::Set< AtomType> s_backbone_atom_types( Begin(), CB.GetIterator());

      // return
      return s_backbone_atom_types;
    }

    //! @brief return StorageVector of backbone atom names
    //! @return StorageVector of backbone atom names
    const storage::Vector< std::string> &AtomTypes::GetBackBoneAtomNames() const
    {
      // initialize vector of strings for storing atom names
      static const storage::Vector< std::string> s_bb_atom_names( Begin(), CB.GetIterator());

      // return
      return s_bb_atom_names;
    }

    //! @brief return StorageVector of side chain Atom Types
    //! @return StorageVector of side chain Atom Types
    const storage::Set< AtomType> &AtomTypes::GetSideChainAtomTypes() const
    {
      // initialize static set of AtomTypes from CB to SG
      static const storage::Set< AtomType> s_side_chain_atom_types( CB.GetIterator(), H.GetIterator());

      // return
      return s_side_chain_atom_types;
    }

    //! @brief return set of AtomTypes composed of CA and CB
    //! @return set of AtomTypes composed of CA and CB
    const storage::Set< AtomType> &AtomTypes::GetCACB() const
    {
      // initialize a set of AtomTypes composed of CA and CB
      static const storage::Set< AtomType> s_ca_cb( CA, CB);

      // return
      return s_ca_cb;
    }

    //! @brief determines and returns the atom type from the provided pdb atom name
    //! @param PDB_ATOM_NAME AtomName for the pdb atom of interest
    //! @return the atom type from the provided pdb atom name
    AtomType AtomTypes::TypeFromPDBAtomName( const std::string &PDB_ATOM_NAME) const
    {
      static std::map< std::string, AtomType> s_map;
      if( s_map.empty())
      {
        for
        (
          const_iterator type_itr( Begin()), type_itr_end( End());
          type_itr != type_itr_end;
          ++type_itr
        )
        {
          s_map[ type_itr->GetName()] = *type_itr;
        }
      }

      auto ite( s_map.find( util::TrimString( PDB_ATOM_NAME)));
      return ite == s_map.end() ? e_Undefined : ite->second;
    }

    //! @brief access to set of possible first side chain atom
    //! @return set of first sied chain atom types, CB and HA2 (for glycin)
    const storage::Set< AtomType> &AtomTypes::GetFirstSidechainAtomTypes() const
    {
      // set of first side chain atom types
      static const storage::Set< AtomType> s_first_sidechain_atom_types( CB, HA2);

      // end
      return s_first_sidechain_atom_types;
    }

    //! @brief access to map of additional atom types for terminal residues and their PDB atom name defined by PDB
    //! @return Map containing the atoms, only found in the terminal amine and carboxylic acid with their PDB atom NAME (e.g. H1 -> 1H)
    const storage::Map< AtomType, std::string> &AtomTypes::GetTerminalExtraAtomTypes() const
    {
      static const storage::Map< AtomType, std::string> s_terminal_atom_type( TerminalExtraAtomTypes());
      return s_terminal_atom_type;
    }

    //! @brief terminal atomtype from atom name
    //! @param ATOM_NAME name of atom (e.g. H1 or 1H, OXT ...)
    //! @return the terminal atom type for that atom name - undefined if there is non
    AtomType AtomTypes::GetTerminalAtomTypeFromName( const std::string &ATOM_NAME) const
    {
      // find atom type by iterating over mapping
      for
      (
        storage::Map< AtomType, std::string>::const_iterator
          itr( GetTerminalExtraAtomTypes().Begin()), itr_end( GetTerminalExtraAtomTypes().End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->second.find( ATOM_NAME) != std::string::npos)
        {
          return itr->first;
        }
      }

      // name can not be found in map
      return TypeFromPDBAtomName( ATOM_NAME);
    }

    //! @brief access to map of additional atom types for terminal residues and their PDB atom name defined by PDB
    //! @return Map containing the atoms, only found in the terminal amine and carboxylic acid with their PDB atom NAME (e.g. H1 -> 1H)
    storage::Map< AtomType, std::string> AtomTypes::TerminalExtraAtomTypes() const
    {
      // set of terminal atom types that only appear at the N and C terminus
      storage::Map< AtomType, std::string> terminal_atom_types;
      terminal_atom_types[  H1] = "1H  ";
      terminal_atom_types[  H2] = "2H  ";
      terminal_atom_types[  H3] = "3H  ";
      terminal_atom_types[ HXT] = " HXT";
      terminal_atom_types[ OXT] = " OXT";

      // end
      return terminal_atom_types;
    }

    //! @brief access to the only instance of AtomTypes
    //! @return reference to only instance of AtomTypes
    const AtomTypes &GetAtomTypes()
    {
      return AtomTypes::GetEnums();
    }

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< biol::AtomTypeData, biol::AtomTypes>;

  } // namespace util
} // namespace bcl
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
#include "biol/bcl_biol_blast_profile.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BlastProfile::BlastProfile() :
      m_Profile( AATypes::s_NumberStandardAATypes, util::GetUndefined< double>()),
      m_Probabilities( AATypes::s_NumberStandardAATypes, util::GetUndefined< double>()),
      m_MatchWeightRelativeToPseudocount( util::GetUndefined< double>()),
      m_Information( util::GetUndefined< double>())
    {
    }
    //! @brief constructor from a profile vector
    //! @param PROFILE blast profile vector
    BlastProfile::BlastProfile
    (
      const linal::Vector< double> &PROFILE,
      const double &MATCH_WEIGHT,
      const double &POSITION_INFORMATION
    ) :
      m_Profile( AATypes::s_NumberStandardAATypes, util::GetUndefined< double>()),
      m_Probabilities( AATypes::s_NumberStandardAATypes, util::GetUndefined< double>()),
      m_MatchWeightRelativeToPseudocount( MATCH_WEIGHT),
      m_Information( POSITION_INFORMATION)
    {
      SetProfile( PROFILE);
    }

    //! @brief virtual copy constructor
    BlastProfile *BlastProfile::Clone() const
    {
      return new BlastProfile( *this);
    }

  ///////////////////////
  // data access - Get //
  ///////////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &BlastProfile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////////////
  // data access - Set //
  ///////////////////////

    //! @brief set BlastProfile to provided PROFILE
    //! @param PROFILE blast profile vector
    void BlastProfile::SetProfile( const linal::Vector< double> &PROFILE)
    {
      // assert the given PROFILE has correct size
      BCL_Assert
      (
        PROFILE.GetSize() == AATypes::s_NumberStandardAATypes,
        "The Blast profile vector passed should have a size of "
        + util::Format().W( 5)( size_t( AATypes::s_NumberStandardAATypes))
        + " but instead it is " + util::Format().W( 5)( PROFILE.GetSize()) + "!!"
      );

      // reset profile
      m_Profile = 0;

      // update profile
      m_Profile = PROFILE;

      // update probabilities
      SetProbabilities();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the conservation, sum of profile times probability
    //! @return the conservation
    double BlastProfile::CalculateConservation() const
    {
      // return the scalar product
      return linal::ScalarProduct( m_Profile, m_Probabilities);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads BlastProfile from ISTREAM
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &BlastProfile::ReadProfile( std::istream &ISTREAM)
    {
      // iterate over profile
      for
      (
        double *ptr( m_Profile.Begin()), *const ptr_end( m_Profile.End());
        ptr != ptr_end;
        ++ptr
      )
      {
        // read value from stream
        ISTREAM >> *ptr;
      }

      // try reading in the probabilities, information per position, and gapless weight
      // some blast implementations (like deltablast) do not provide this information, so set probabilities up from
      // blast profile if the probabilities are not availability.
      std::string rest_of_line;
      std::getline( ISTREAM, rest_of_line);
      storage::Vector< double> remaining_numbers( util::SplitStringToNumerical< double>( rest_of_line));

      // four possibilities: blast probabilities given, information gain and match weight given, both, or neither
      if( ( remaining_numbers.GetSize() % 20) == size_t( 2))
      {
        // information and alignment weight were given, read them in
        m_Information = remaining_numbers( remaining_numbers.GetSize() - 2);
        m_MatchWeightRelativeToPseudocount = remaining_numbers( remaining_numbers.GetSize() - 1);
      }

      // update probabilities
      SetProbabilities();

      // return
      return ISTREAM;
    }

    //! @brief writes BlastProfile to OSTREAM
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &BlastProfile::WriteProfile( std::ostream &OSTREAM) const
    {
      // iterate over profile
      util::Format format;
      format.W( 7).FFP( 2);
      for
      (
        const double *ptr( m_Profile.Begin()), *const ptr_end( m_Profile.End());
        ptr != ptr_end;
        ++ptr
      )
      {
        // write value to stream
        OSTREAM << format( *ptr);
      }
      OSTREAM << format( m_Information) << format( m_MatchWeightRelativeToPseudocount);
      // return
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &BlastProfile::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Profile, ISTREAM);
      io::Serialize::Read( m_Probabilities, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &BlastProfile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // read members
      io::Serialize::Write( m_Profile, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Probabilities, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

    //! @brief set BlastProfile probabilities
    void BlastProfile::SetProbabilities()
    {
      // calculate probabilities by taking 10 ^ ( m_profile/10)
      m_Probabilities = double( 10) ^ ( m_Profile / double( 10));

      // normalize (typically by 20)
      m_Probabilities /= m_Probabilities.Sum();
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_blast_profile_handler.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BlastProfileHandler::BlastProfileHandler()
    {
    }

    //! @brief virtual copy constructor
    BlastProfileHandler *BlastProfileHandler::Clone() const
    {
      return new BlastProfileHandler( *this);
    }

    //! @brief virtual destructor
    BlastProfileHandler::~BlastProfileHandler()
    {
    }

  ///////////////////////
  // data access - Get //
  ///////////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &BlastProfileHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read BLAST profile from ISTREAM for given amino acid
    //! @param ISTREAM input stream
    //! @param AMINO_ACID AABase into which blast profile is going to read
    //! @return std::istream which was read from
    std::istream &BlastProfileHandler::ReadProfileForAA( std::istream &ISTREAM, AABase &AMINO_ACID)
    {
      BCL_Assert( ISTREAM.good(), "Truncated blast file; aa: " + AMINO_ACID.GetIdentification());
      // read seqid and amino acid one letter code
      int seqid;
      std::string one_letter_code;
      ISTREAM >> seqid >> one_letter_code;

      // get current aa type from read one letter code
      const AAType current_aa_type( GetAATypes().AATypeFromOneLetterCode( one_letter_code[ 0]));

      // check matching amino acid seqids
      BCL_Assert
      (
        AMINO_ACID.GetSeqID() == seqid,
        "mismatch in seqids! sequence: " + AMINO_ACID.GetIdentification() +
        " vs. from blast: " + util::Format()( seqid) + " " + current_aa_type.GetName()
      );

      // if types of residue from sequence and from file do not match
      // and amino acid in the sequence is arbitrary/undefined
      if( current_aa_type != AMINO_ACID.GetType())
      {
        if
        (
          !AMINO_ACID.GetType().IsDefined()
          || AMINO_ACID.GetType() == GetAATypes().XXX
          || AMINO_ACID.GetType() == GetAATypes().UNK
        )
        {
          // BLAST has a valid AA type but the sequence does not.  Warn the user that we're about to trust BLAST over the
          // PDB file
          BCL_MessageCrt
          (
            " changing the type of the aa " + AMINO_ACID.GetIdentification() + " to " + current_aa_type.GetName()
            + " to match blast profile"
          );
          // change aa type
          AMINO_ACID.SetData
          (
            util::ShPtr< AAData>
            (
              new AAData
              (
                current_aa_type,
                AMINO_ACID.GetSeqID(),
                AMINO_ACID.GetPdbID(),
                AMINO_ACID.GetPdbICode(),
                AMINO_ACID.GetChainID()
               )
            )
          );
        }
        else if
        (
          !current_aa_type.IsDefined()
          || current_aa_type == GetAATypes().UNK
        )
        {
          // blast profile had an undefined type; make sure that this corresponds to an unnatural AA in the sequence
          // or that the convert to natural aa type flag is set
          BCL_Assert
          (
            !AMINO_ACID.GetType()->IsNaturalAminoAcid() || pdb::Factory::GetFlagConvertToNaturalAAType()->GetFlag(),
            "mismatch in amino acid types ! sequence: chain: " + util::Format()( AMINO_ACID.GetChainID()) + " " +
            AMINO_ACID.GetIdentification() +
            " vs. from blast: " + util::Format()( seqid) + " " + current_aa_type.GetName()
          );
        }
        else if( current_aa_type == GetAATypes().XXX)
        {
          // Indicates a low complexity or tandem repeat region when PSSM is generated by a program that hard-masks such
          // regions (e.g. tantan -x X used with psiblast)
        }
        else if( current_aa_type->GetParentType() != AMINO_ACID.GetType()->GetParentType())
        {
          // bcl and blast profile type different -> real problem
          BCL_Exit
          (
            "mismatch in amino acid types ! sequence: chain: " + util::Format()( AMINO_ACID.GetChainID()) + " " +
            AMINO_ACID.GetIdentification() +
            " vs. from blast: " + util::Format()( seqid) + " " + current_aa_type.GetName(),
            -1
          );
        }
      }

      // set the profile to read vector
      AMINO_ACID.ReadBlastProfile( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write BLAST profile to OSTREAM for given amino acid
    //! @param OSTREAM output stream
    //! @param AMINO_ACID AABase f which blast profile is going to read
    //! @return std::ostream which was written to
    std::ostream &BlastProfileHandler::WriteProfileForAA( std::ostream &OSTREAM, const AABase &AMINO_ACID)
    {
      // output the seq id and one letter code
      OSTREAM << util::Format().W( 5)( AMINO_ACID.GetSeqID()) << ' ' << AMINO_ACID.GetType()->GetOneLetterCode() << ' ';

      // output blast profile
      AMINO_ACID.GetBlastProfile().WriteProfile( OSTREAM);

      // write line break
      OSTREAM << '\n';

      // return
      return OSTREAM;
    }

    //! @brief read BLAST profile from ISTREAM for given AASequence
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which blast profile is going to read
    //! @return std::istream which was read from
    std::istream &BlastProfileHandler::ReadProfileForAASequence
    (
      std::istream &ISTREAM,
      AASequence &AA_SEQUENCE
    )
    {
      // local variables
      std::string line;

      // read line by line until beginning of blast profile is found
      while( std::getline( ISTREAM, line) && ( !line.length() || line[ line.length() - 1] != 'V'))
      {
      }

      // iterate over amino acids
      for
      (
        AASequence::iterator aa_itr( AA_SEQUENCE.Begin()), aa_itr_end( AA_SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // read the Blast profile for this amino acid from ISTREAM
        ReadProfileForAA( ISTREAM, **aa_itr);
      }

      // end
      return ISTREAM;
    }

    //! @brief write BLAST profile to OSTREAM for given AASequence
    //! @param OSTREAM output stream
    //! @param AA_SEQUENCE AASequence from which blast profile is going to read
    //! @return std::ostream which was written to
    std::ostream &BlastProfileHandler::WriteProfileForAASequence
    (
      std::ostream &OSTREAM,
      const AASequence &AA_SEQUENCE
    )
    {
      // write header to OSTREAM
      OSTREAM << "\n# This file was rebuilt for "
              << AA_SEQUENCE.GetFastaHeader()
              << " using AASequence::WriteBlastProfile\n        ";

      // iterate over 20 standard amino acid type to form the first line of the blast profile
      for
      (
        AATypes::const_iterator aa_itr( GetAATypes().Begin()),
          aa_itr_end( GetAATypes().GetEnumIteratorFromIndex( AATypes::s_NumberStandardAATypes));
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // output the one letter code for this amino acid in the first line
        OSTREAM << "  " << ( *aa_itr)->GetOneLetterCode();
      }

      // output a line break
      OSTREAM << '\n';

      // iterate over amino acids in the sequence
      for
      (
        AASequence::const_iterator aa_itr( AA_SEQUENCE.Begin()), aa_itr_end( AA_SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // output the blast profile for this amino acid to OSTREAM
        WriteProfileForAA( OSTREAM, **aa_itr);
      }

      // end
      return OSTREAM;

    }

    //! @brief read BLAST profile from ISTREAM for given ProteinModel
    //! @param PROTEIN_MODEL ProteinModel into which blast profile is going to read
    //! @param PREFIX prefix of the protein, including path ( usually pdb id)
    //! @param EXTENSION the desired extension; if blank, .ascii6 will be preferred, but .ascii will also be accepted
    //! @return whether reading was successful
    bool BlastProfileHandler::TryReadProfileForProteinModel
    (
      assemble::ProteinModel &PROTEIN_MODEL,
      const std::string &PREFIX,
      const std::string &EXTENSION
    )
    {
      // initialize success booelean to true
      bool success( true);

      // valid strings to use if the chain id is blank
      const storage::Vector< std::string> s_blank_chain_id_names( storage::Vector< std::string>::Create( " ", "_", ""));

      // create primary and alternative extensions
      const std::string pref_extension( EXTENSION.empty() ? ".ascii6" : EXTENSION);
      const std::string alt_extension( EXTENSION.empty() ? ".ascii" : EXTENSION);

      //iterate over all chain and read blast profile from file generated from PATH, SEQ_TAG and ChainID
      for
      (
        util::ShPtrVector< assemble::Chain>::iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // test whether this chain already has a valid blast profile, if so, skip it
        if
        (
          ( *chain_itr)->GetSequence()->GetSize() == size_t( 0)
          || ( *chain_itr)->GetSequence()->GetFirstAA()->GetBlastProfilePtr().IsDefined()
        )
        {
          BCL_MessageDbg( "Skipping reading blast profile for " + PREFIX + "; already read. ");
          continue;
        }

        const std::string chain_id( size_t( 1), ( *chain_itr)->GetChainID());

        // store the valid path in this string
        std::string blast_path;

        // handle blank chains
        if( io::DirectoryEntry( PREFIX + chain_id + pref_extension).DoesExist())
        {
          blast_path = PREFIX + chain_id + pref_extension;
        }
        else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + chain_id + alt_extension).DoesExist())
        {
          blast_path = PREFIX + chain_id + alt_extension;
        }
        else if( io::DirectoryEntry( PREFIX + pref_extension).DoesExist())
        {
          blast_path = PREFIX + pref_extension;
        }
        else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + alt_extension).DoesExist())
        {
          blast_path = PREFIX + alt_extension;
        }
        else if( chain_id[ 0] == ' ') // handle blank / unknown chains
        {
          if( io::DirectoryEntry( PREFIX + "_" + pref_extension).DoesExist())
          {
            blast_path = PREFIX + "_" + pref_extension;
          }
          else if( io::DirectoryEntry( PREFIX + pref_extension).DoesExist())
          {
            blast_path = PREFIX + pref_extension;
          }
          else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + "_" + alt_extension).DoesExist())
          {
            blast_path = PREFIX + "_" + alt_extension;
          }
          else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + alt_extension).DoesExist())
          {
            blast_path = PREFIX + alt_extension;
          }
        }

        if( blast_path.empty()) // was any valid blast path found?
        {
          // no valid blast path found
          BCL_MessageVrb
          (
            "Failed to find blast profile for chain '" + chain_id + "' with prefix: " + PREFIX
          );
          success = false;
        }
        else
        {
          io::IFStream input;
          success = io::File::TryOpenIFStream( input, blast_path) && success;
          if( success)
          {
            BCL_MessageDbg
            (
              "Reading Blastprofile for Chain '" + chain_id + "' from " + blast_path
            );

            // read blast profile for this chain
            ReadProfileForAASequence( input, *( *chain_itr)->GetSequence());
            // clear stream
            io::File::CloseClearFStream( input);
          }
        }
      }

      //return success
      return success;
    }

    //! @brief read all the blast profiles for a given sequence into a series of vectors.
    //!        This is the only method of reading that does not change the underlying data structures and is likewise
    //!        safe inside threaded code that may be operating on the same sequence
    //! @param SEQUENCE the sequence of interest; only used to check that the blast profile has the same sequence
    //! @param PREFIX prefix of the protein, including path ( usually pdb id)
    //! @param EXTENSION the desired extension; if blank, .ascii6 will be preferred, but .ascii will also be accepted
    storage::Vector< BlastProfile> BlastProfileHandler::ReadProfilesForConstAASequence
    (
      const AASequence &SEQUENCE,
      const std::string &PREFIX,
      const std::string &EXTENSION
    )
    {
      storage::Vector< BlastProfile> blast_profiles;

      // create primary and alternative extensions
      const std::string pref_extension( EXTENSION.empty() ? ".ascii6" : EXTENSION);
      const std::string alt_extension( EXTENSION.empty() ? ".ascii" : EXTENSION);
      io::IFStream input;

      // test whether this chain is empty if so, skip it
      if( SEQUENCE.GetSize() == size_t( 0))
      {
        BCL_MessageDbg( "Skipping reading blast profile for empty sequence with prefix " + PREFIX + "");
        return blast_profiles;
      }

      const std::string chain_id( size_t( 1), SEQUENCE.GetChainID());

      // store the valid path in this string
      std::string blast_path;

      // handle blank chains
      if( io::DirectoryEntry( PREFIX + chain_id + pref_extension).DoesExist())
      {
        blast_path = PREFIX + chain_id + pref_extension;
      }
      else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + chain_id + alt_extension).DoesExist())
      {
        blast_path = PREFIX + chain_id + alt_extension;
      }
      else if( io::DirectoryEntry( PREFIX + pref_extension).DoesExist())
      {
        blast_path = PREFIX + pref_extension;
      }
      else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + alt_extension).DoesExist())
      {
        blast_path = PREFIX + alt_extension;
      }
      else if( chain_id[ 0] == ' ') // handle blank / unknown chains
      {
        if( io::DirectoryEntry( PREFIX + "_" + pref_extension).DoesExist())
        {
          blast_path = PREFIX + "_" + pref_extension;
        }
        else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + "_" + alt_extension).DoesExist())
        {
          blast_path = PREFIX + "_" + alt_extension;
        }
      }

      if( blast_path.empty()) // was any valid blast path found?
      {
        // no valid blast path found
        BCL_MessageCrt
        (
          "Failed to find blast profile for chain '" + chain_id + "' with prefix: " + PREFIX + " and suffix: "
          + pref_extension + ( pref_extension == alt_extension ? "" : " or " + alt_extension)
        );
        return blast_profiles;
      }

      io::File::MustOpenIFStream( input, blast_path);
      BCL_MessageDbg( "Reading Blastprofile for Chain '" + chain_id + "' from " + blast_path);

      // read blast profile for this chain
      ReadProfilesForConstAASequence( SEQUENCE, blast_profiles, input);
      // clear stream
      io::File::CloseClearFStream( input);

      //return success
      return blast_profiles;
    }

    //! @brief read all the blast profiles for a given sequence into a series of blast profiles
    //!        This is the only method of reading that does not change the underlying data structures and is likewise
    //!        safe inside threaded code that may be operating on the same sequence
    //! @param SEQUENCE the sequence of interest; only used to check that the blast profile has the same sequence
    //! @param PROFILES reference to a vector that will be used to store the profiles
    //! @param ISTREAM input stream
    void BlastProfileHandler::ReadProfilesForConstAASequence
    (
      const AASequence &SEQUENCE,
      storage::Vector< BlastProfile> &PROFILES,
      std::istream &ISTREAM
    )
    {
      PROFILES.Reset();
      PROFILES.AllocateMemory( SEQUENCE.GetSize());

      // read line by line until beginning of blast profile is found
      for( std::string line; std::getline( ISTREAM, line) && ( !line.length() || line[ line.length() - 1] != 'V');)
      {
      }

      // iterate over amino acids
      for
      (
        AASequence::const_iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        const AABase &actual_aa( **aa_itr);
        util::ShPtr< AAData> copy_aadata( new AAData( actual_aa.GetType(), actual_aa.GetSeqID()));
        // create an aa
        AA copy_aa( copy_aadata);
        // copy the aa
        // read the Blast profile for this amino acid from ISTREAM
        ReadProfileForAA( ISTREAM, copy_aa);
        PROFILES.PushBack( copy_aa.GetBlastProfile());
      }
    }

    //! @brief write BLAST profile to OSTREAM for given ProteinModel
    //! @param PROTEIN_MODEL ProteinModel from which blast profile is going to read
    //! @param PREFIX prefix of the protein ( usually pdb id)
    //! @param PATH the directory path where blast files will be written
    //! @return whether writing was successful
    bool BlastProfileHandler::WriteProfileForProteinModel
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const std::string &PREFIX,
      const std::string &PATH
    )
    {
      // iterate over chains in the model
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // create the output file for this chain's blast profile
        std::string output_filename
        (
          PATH + PATH_SEPARATOR + PREFIX + ( *chain_itr)->GetChainID() + ".ascii"
        );

        // create ofstream and check it
        io::OFStream write( output_filename.c_str());
        BCL_Assert( write, "The file cannot be opened for output " + output_filename);

        // call write function of the chain
        WriteProfileForAASequence( write, *( *chain_itr)->GetSequence());

        // reset ofstream
        io::File::CloseClearFStream( write);

        // return
        return true;
      }

      // return
      return true;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &BlastProfileHandler::Read( std::istream &ISTREAM)
    {
      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &BlastProfileHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_chi_angle.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
  //////////
  // data //
  //////////

    //! @brief conversion to a string from a Chi
    //! @param CHI the chi angle to get a string for
    //! @return a string representing that chi
    const std::string &ChiAngle::GetChiName( const ChiAngle::Chi &CHI)
    {
      static const std::string s_descriptors[] =
      {
        "e_One",
        "e_Two",
        "e_Three",
        "e_Four",
        "e_Five",
        "e_Undefined",
        GetStaticClassName< ChiAngle>()
      };
      return s_descriptors[ size_t( CHI)];
    }

      //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ChiAngle::s_Instance
    (
      GetObjectInstances().AddInstance( new ChiAngle())
    );

    //! @brief gives the maximum magnitude that a chi angle could have
    //! @param UNIT the unit the max angle should be given in
    //! @return double which is the maximum magnitude that a chi angle could have
    double ChiAngle::GetMaxAngle( const math::Angle::Unit &UNIT)
    {
      if( UNIT == math::Angle::e_Degree)
      {
        return 180;
      }
      if( UNIT == math::Angle::e_Radian)
      {
        return math::g_Pi;
      }

      return util::GetUndefinedDouble();
    }

    //! @brief gives the total number of angular units in a circle
    //! @param UNIT the unit the angular units are in
    //! @return double the total number of angular units in a circle
    double ChiAngle::GetCircle( const math::Angle::Unit &UNIT)
    {
      return GetMaxAngle( UNIT) * 2.0;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ChiAngle::ChiAngle() :
      m_Chi( ChiAngle::e_Undefined),
      m_Angle( util::GetUndefinedDouble()),
      m_Unit( math::Angle::e_Undefined)
    {
    }

    //! @brief constructor from member variable parameters
    //! @param CHI the chi this will correspond to
    ChiAngle::ChiAngle( const Chi &CHI) :
      m_Chi( CHI),
      m_Angle( util::GetUndefinedDouble()),
      m_Unit( math::Angle::e_Undefined)
    {
    }

    //! @brief constructor from member variable parameters
    //! @param CHI the chi this will correspond to
    //! @param ANGLE the angle of this chi
    //! @param UNIT the unit type the angle is measured in
    ChiAngle::ChiAngle( const Chi &CHI, const double ANGLE, const math::Angle::Unit &UNIT) :
      m_Chi( CHI),
      m_Angle( ANGLE),
      m_Unit( UNIT)
    {
      BCL_Assert
      (
        math::Absolute( GetAngle( m_Unit) <= GetMaxAngle( m_Unit)) || !util::IsDefined( ANGLE),
        "chi angle in units of " + m_Unit.GetString() +
        " should not be larger in magnitude than " + util::Format()( GetMaxAngle( m_Unit))
      );
    }

    //! @brief Clone function
    //! @return pointer to new ChiAngle
    ChiAngle *ChiAngle::Clone() const
    {
      return new ChiAngle( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ChiAngle::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the chi this corresponds to
    //! @return Chi which is the chi this corresponds to
    const ChiAngle::Chi &ChiAngle::GetChi() const
    {
      return m_Chi;
    }

    //! @brief gives the angular value of this angle
    //! @param ANGLE_UNIT the unit the angle should be given in
    //! @return double which is the angular value of this angle in the units desired according ANGLE_UNIT
    double ChiAngle::GetAngle( const math::Angle::Unit &ANGLE_UNIT) const
    {
      // true if angle stored as degrees and degrees are desired
      if( m_Unit == math::Angle::e_Degree && ANGLE_UNIT == math::Angle::e_Degree)
      {
        // return the angle
        return m_Angle;
      }
      // true if angle stored as radians and degrees are desired
      if( m_Unit == math::Angle::e_Radian && ANGLE_UNIT == math::Angle::e_Degree)
      {
        // convert the radians to degrees and return value
        return math::Angle::Degree( m_Angle);
      }
      // true if angle stored as degrees and radians are desired
      if( m_Unit == math::Angle::e_Degree && ANGLE_UNIT == math::Angle::e_Radian)
      {
        // convert the degrees to radian and return value
        return math::Angle::Radian( m_Angle);
      }
      // true if angle stored as radians and radians are desired
      if( m_Unit == math::Angle::e_Radian && ANGLE_UNIT == math::Angle::e_Radian)
      {
        // return the angle
        return m_Angle;
      }

      return util::GetUndefinedDouble();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates the absolute angular difference between this chi angle and a provided chi angle
    //! @param CHI_ANGLE the chi angle whose difference will be calculated from this
    //! @param ANGLE_UNIT the unit the angular difference should be given in
    //! @return double which is the absolute angular difference between this and the given ChiAngle in desired units
    double ChiAngle::CalculateAngleDifference( const ChiAngle &CHI_ANGLE, const math::Angle::Unit &ANGLE_UNIT) const
    {
      // true if the chi angles are not for the same chi or either of the chi angles are not defined
      if
      (
        CHI_ANGLE.GetChi() != GetChi() || !util::IsDefined( CHI_ANGLE.GetAngle( ANGLE_UNIT))
        || !util::IsDefined( GetAngle( ANGLE_UNIT))
      )
      {
        return util::GetUndefinedDouble();
      }

      double absolute_difference( math::Absolute( CHI_ANGLE.GetAngle( ANGLE_UNIT) - GetAngle( ANGLE_UNIT)));

      // need to take into account that if one angle is -175 and the other is
      // 175 the difference is 10 degrees, not 350
      // so if the difference is >= than 180 then subtract it from 360
      // in order to get the smaller portion of the circle between the angles
      if( absolute_difference > GetMaxAngle( ANGLE_UNIT))
      {
        absolute_difference = GetCircle( ANGLE_UNIT) - absolute_difference;
      }

      return absolute_difference;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ChiAngle::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Chi, ISTREAM);
      io::Serialize::Read( m_Angle, ISTREAM);
      io::Serialize::Read( m_Unit, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ChiAngle::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Chi, OSTREAM, INDENT);
      io::Serialize::Write( m_Angle, OSTREAM, INDENT);
      io::Serialize::Write( m_Unit, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief reads stream formatted in way easier for user to create
    //!        format is space separated
    //!        <chi enum> <angle value> <angle unit>
    //!        An example is
    //!        e_Two 90 degree
    //!        Another example is
    //!        e_Three 1.1 radian
    //! @return istream the ChiAngle was read from
    std::istream &ChiAngle::ReadSimple( std::istream &ISTREAM)
    {
      std::string chi;
      std::getline( ISTREAM, chi);
      BCL_MessageDbg( "chi line is |" + chi + "|");
      if( !chi.empty())
      {
        std::stringstream read( chi);
        io::Serialize::Read( m_Chi, read);
        io::Serialize::Read( m_Angle, read);
        io::Serialize::Read( m_Unit, read);
      }

      return ISTREAM;
    }

    //! @brief gives description of this in format as read by ReadSimple
    //! @return std::string gives description of this in ReadSimple format
    std::string ChiAngle::WriteSimple( const math::Angle::Unit &ANGLE_UNIT) const
    {
      return m_Chi.GetString() + " " +
      util::Format().W( 8).FFP( 3)( GetAngle( ANGLE_UNIT)) + " " + math::Angle::UnitEnum( ANGLE_UNIT).GetString();
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace biol
} // namespace bcl
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

// This code is adapted from Pteros; see http://pteros.sourceforge.net/dssp_8h_source.html
// Pteros' copyright, applies solely to this file:
// Portions Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
// Boost Software License - Version 1.0 - August 17th, 2003
//
// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:
//
// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "biol/bcl_biol_dssp.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_back_bone_completer.h"
#include "biol/bcl_biol_atom.h"
#include "math/bcl_math_mutate_result.h"
#include "sspred/bcl_sspred_pdb.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief compare bridge by chainid and seqid of residues in bridge
    //! @param BRIDGE bridge to comparet his to
    //! @return true if this bridge comes before the argument BRIDGE in sequence
    bool DSSP::Bridge::operator<( const Bridge &BRIDGE) const
    {
      return m_ChainI < BRIDGE.m_ChainI || ( m_ChainI == BRIDGE.m_ChainI && m_I.front()->GetSeqID() < BRIDGE.m_I.front()->GetSeqID());
    }

    //! @brief constructor from members
    DSSP::BridgePartner::BridgePartner() :
      m_Partner(),
      m_Ladder(),
      m_Parallel()
    {
    }

    //! @brief constructor from members
    DSSP::BridgePartner::BridgePartner( const util::SiPtr< const AABase> &PARTNER, const size_t LADDER, const bool PARALLEL) :
      m_Partner( PARTNER),
      m_Ladder( LADDER),
      m_Parallel( PARALLEL)
    {
    }

  //////////
  // data //
  //////////

    //! @brief default max HBond energy for an HBond to be considered
    //! @return default max energy
    const double &DSSP::GetDefaultHBondMaxEnergy()
    {
      static const double s_max_hydrogen_bond_energy( -0.45);
      return s_max_hydrogen_bond_energy;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param MAX_HYDROGEN_BOND_ENERGY maximal hydrogen bond energy to be considered an H-bond
    DSSP::DSSP( const double MAX_HYDROGEN_BOND_ENERGY) :
      m_MaxHydrogenBondEnergy( MAX_HYDROGEN_BOND_ENERGY),
      m_NrOfHydrogenBondsInParallelBridges( 0),
      m_NrOfHydrogenBondsInAntiparallelBridges( 0)
    {
      Reset();
    }

    //! @brief Clone function
    //! @return pointer to new DSSP
    DSSP *DSSP::Clone() const
    {
      return new DSSP( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DSSP::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &DSSP::GetScheme() const
    {
      static const std::string s_scheme( "dssp");
      return s_scheme;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get various statistics
    void DSSP::GetStatistics
    (
      size_t &NR_OF_RESIDUES,
      size_t &NR_OF_CHAINS,
      size_t &NR_OF_SS_BRIDGES,
      size_t &NR_OF_INTRA_CHAIN_SS_BRIDGES,
      size_t &NR_OF_HYDROGEN_BONDS,
      size_t NR_OF_HYDROGEN_BONDS_PER_DISTANCE[ s_HistogramSize],
      size_t NR_RES_PER_ALPHA_HELIX_HISTOGRAM[ s_HistogramSize]
    ) const
    {
      const double max_hydrogen_bond_energy( -0.5);
      NR_OF_RESIDUES = m_SecondaryStructure.size();
      NR_OF_CHAINS   = 0;
      NR_OF_SS_BRIDGES = 0;
      NR_OF_INTRA_CHAIN_SS_BRIDGES = 0;
      NR_OF_HYDROGEN_BONDS = 0;
      std::fill_n( NR_OF_HYDROGEN_BONDS_PER_DISTANCE, 11, 0.0);

      for
      (
        HydrogenBondContainerType::const_iterator
          itr( m_HydrogenBondAcceptor.Begin()), itr_end( m_HydrogenBondAcceptor.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->second.First().Second() < max_hydrogen_bond_energy)
        {
          ++NR_OF_HYDROGEN_BONDS;
          const int seq_sep( itr->first->GetSeqID() - itr->second.First().First()->GetSeqID());
          if( seq_sep >= -5 && seq_sep <= 5 && itr->first->GetChainID() == itr->second.First().First()->GetChainID())
          {
            ++NR_OF_HYDROGEN_BONDS_PER_DISTANCE[ seq_sep + 5];
          }
        }
      }

      // iterate through chains
      std::fill_n( NR_RES_PER_ALPHA_HELIX_HISTOGRAM, s_HistogramSize, 0);
      for
      (
        assemble::ProteinModel::const_iterator chain_itr( m_Model->GetChains().Begin()), chain_itr_end( m_Model->GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        size_t helix_length( 0);
        const util::SiPtrVector< const AABase> res( ( *chain_itr)->GetAminoAcids());
        for
        (
          util::SiPtrVector< const AABase>::const_iterator aa_itr( res.Begin()), aa_itr_end( res.End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          if( m_SecondaryStructure.find( *aa_itr)->second == e_Alphahelix)
          {
            ++helix_length;
          }
          else if( helix_length > 0)
          {
            helix_length = std::min( helix_length, size_t( s_HistogramSize));
            NR_RES_PER_ALPHA_HELIX_HISTOGRAM[ helix_length - 1] += 1;
            helix_length = 0;
          }
        }
      }
    }

    //! @brief reset all members
    void DSSP::Reset() const
    {
      m_Model = util::ShPtr< assemble::ProteinModel>();

      m_HydrogenBondDonor.Reset();
      m_HydrogenBondAcceptor.Reset();

      m_BridgePartner.Reset();

      m_Sheet.Reset();

      m_SecondaryStructure.clear();

      for( size_t stride( 3); stride <= 5; ++stride)
      {
        m_HelixFlag[ stride - 3].clear();
      }

      m_Bend.Reset();

      m_Alpha.clear();

      m_Kappa.clear();
      m_Phi.clear();
      m_Psi.clear();
      m_TCO.clear();

      m_NrOfHydrogenBondsInParallelBridges = 0;
      m_NrOfHydrogenBondsInAntiparallelBridges = 0;
      std::fill_n( m_ParallelBridgesPerLadderHistogram, s_HistogramSize, 0);
      std::fill_n( m_AntiparallelBridgesPerLadderHistogram, s_HistogramSize, 0);
      std::fill_n( m_LaddersPerSheetHistogram, s_HistogramSize, 0);
    }

    //! @brief Set PDB SS Predictions using the DSSP algorithm
    //! @param PROTEIN_MODEL the protein model with sses (can be just one single loop)
    void DSSP::SetPDBSSPred( assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // construct a new model with the dssp algorithm
      util::ShPtr< assemble::ProteinModel> sp_dssp_model( operator()( PROTEIN_MODEL).GetArgument());

      // iterate through the chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          dssp_chain_itr( sp_dssp_model->GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr, ++dssp_chain_itr
      )
      {
        // get all the sses from the dssp-fixed chain
        util::SiPtrVector< const assemble::SSE> dssp_sses( ( *dssp_chain_itr)->GetSSEs());

        // make a map from AASeqID to dssp-related SSType
        storage::Map< int, SSType> seq_id_to_dssp_type;
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator
            sse_itr( dssp_sses.Begin()), sse_itr_end( dssp_sses.End());
          sse_itr != sse_itr_end;
          ++sse_itr
        )
        {
          const assemble::SSE &dssp_sse( **sse_itr);
          // get the type of this sse
          SSType ss_type( dssp_sse.GetType());
          // get start and finish

          for
          (
            int seq_id( dssp_sse.GetFirstAA()->GetSeqID()), end_seq_id( dssp_sse.GetLastAA()->GetSeqID() + 1);
            seq_id < end_seq_id;
            ++seq_id
          )
          {
            seq_id_to_dssp_type[ seq_id] = ss_type;
          }
        }

        // iterate through the sequence
        for
        (
          AASequence::const_iterator
            aa_itr( ( *chain_itr)->GetSequence()->Begin()),
            aa_itr_end( ( *chain_itr)->GetSequence()->End());
          aa_itr != aa_itr_end; ++aa_itr
        )
        {
          // create a pointer to the data
          util::ShPtr< AABase> sp_aa( *aa_itr);
          const util::SiPtr< const sspred::MethodInterface> sp_aa_ss_pdb( sp_aa->GetSSPrediction( sspred::GetMethods().e_PDB));
          SSType existing_ss_type( GetSSTypes().COIL), ss_type( GetSSTypes().COIL);
          EnvironmentType existing_environment( GetEnvironmentTypes().e_Solution);
          if( !sp_aa_ss_pdb.IsDefined())
          {
            BCL_MessageCrt
            (
              "Warning: AA with id: " + sp_aa->GetIdentification() + " had no secondary structure from the pdb file, "
              "setting it to type coil"
            );
          }
          else
          {
            existing_ss_type = sp_aa_ss_pdb->GetOneStateSSPrediction();
            existing_environment = sp_aa_ss_pdb->GetOneStateTMPrediction();
          }

          // look for this seq id in the map
          storage::Map< int, SSType>::const_iterator itr_type_map( seq_id_to_dssp_type.Find( sp_aa->GetSeqID()));
          if( itr_type_map == seq_id_to_dssp_type.End())
          {
            // no type, set it to coil
            ss_type = GetSSTypes().COIL;
          }
          else
          {
            ss_type = itr_type_map->second;
          }
          if( existing_ss_type != ss_type)
          {
            // update the secondary structure prediction for the AA
            sp_aa->SetSSPrediction( sspred::GetMethods().e_PDB, sspred::PDB( ss_type, existing_environment));
          }
        }
      }
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that uses the coordinates of the secondary structure elements within them to generate a new
    //!        new secondary structure assignment using DSSP
    //! @param PROTEIN_MODEL the protein model with sses (can be just one single loop)
    //! @return a mutate result with a protein model that has new secondary structure elements based on dssp
    math::MutateResult< assemble::ProteinModel> DSSP::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      Reset();
      AABackBoneCompleter completer( true, false, false);
      m_Model = completer.CompleteProteinModel( PROTEIN_MODEL);

      CalculateHydrogenBondEnergies();
      CalculateBetaSheets();
      CalculateAlphaHelices();

      util::ShPtrVector< assemble::Chain> new_chains;

      for
      (
        assemble::ProteinModel::const_iterator chain_itr( m_Model->GetChains().Begin()), chain_itr_end( m_Model->GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        util::ShPtr< assemble::Chain> sp_chain( new assemble::Chain( ( *chain_itr)->GetSequence()));

        const util::SiPtrVector< const AABase> res( ( *chain_itr)->GetAminoAcids());

        util::ShPtrVector< AABase> new_aas;
        SSType last_ss_type( GetSSTypes().COIL);

        int last_seq_id( 0);

        for
        (
          util::SiPtrVector< const AABase>::const_iterator aa_itr( res.Begin()), aa_itr_end( res.End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          const SecondaryStructureType ss_type( m_SecondaryStructure.find( *aa_itr)->second);
          SSType new_ss_type( GetSSTypes().COIL);
          switch( ss_type)
          {
            case e_Loop:       new_ss_type = GetSSTypes().COIL; break;
            case e_Alphahelix: new_ss_type = GetSSTypes().HELIX; break;
            case e_Helix_3:    new_ss_type = GetSSTypes().e_HelixRight310; break;
            case e_Helix_5:    new_ss_type = GetSSTypes().e_HelixRightPi; break;
            case e_Strand:     new_ss_type = GetSSTypes().STRAND; break;
            default:           new_ss_type = GetSSTypes().COIL; break;
          }
          if
          (
               ( new_ss_type == last_ss_type && last_seq_id + 1 == ( *aa_itr)->GetSeqID())
            || new_aas.IsEmpty()
          )
          {
            new_aas.PushBack( util::ShPtr< AABase>( ( *aa_itr)->Clone()));
          }
          else
          {
            util::ShPtr< assemble::SSE> sp_new_sse( new assemble::SSE( AASequence( new_aas, sp_chain->GetChainID()), last_ss_type));
            sp_chain->Insert( sp_new_sse);
            new_aas.Reset();
            new_aas.PushBack( util::ShPtr< AABase>( ( *aa_itr)->Clone()));
          }

          last_ss_type = new_ss_type;
          last_seq_id = ( *aa_itr)->GetSeqID();
        } // amino acids

        if( !new_aas.IsEmpty())
        {
          util::ShPtr< assemble::SSE> sp_new_sse( new assemble::SSE( AASequence( new_aas, sp_chain->GetChainID()), last_ss_type));
          sp_chain->Insert( sp_new_sse);
        }

        new_chains.PushBack( sp_chain);
      } // chains

      util::ShPtr< assemble::ProteinModel> sp_new_model( new assemble::ProteinModel( new_chains));

      return math::MutateResult< assemble::ProteinModel>( sp_new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to file in standard DSSP format
    //! @param OSTREAM stream to write to
    //! @return the ostream written to
    std::ostream &DSSP::WriteToFile( std::ostream &OSTREAM) const
    {
      OSTREAM << "==== Secondary Structure Definition by the program DSSP, BCL version ==== \n";

      size_t nr_of_residues;
      size_t nr_of_chains;
      size_t nr_of_ss_bridges;
      size_t nr_of_intra_chain_ss_bridges;
      size_t nr_of_hydrogen_bonds;
      size_t nr_of_hydrogen_bonds_per_distance[ 11] = {};

      size_t nr_res_per_alpha_helix_histogram[ s_HistogramSize];
      GetStatistics
      (
        nr_of_residues,
        nr_of_chains,
        nr_of_ss_bridges,
        nr_of_intra_chain_ss_bridges,
        nr_of_hydrogen_bonds,
        nr_of_hydrogen_bonds_per_distance,
        nr_res_per_alpha_helix_histogram
      );

      OSTREAM << "REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637\n";

      const util::Format count_format( util::Format().W( 3).Fill( ' ').R());
      const util::Format per100_format( util::Format().W( 5).Fill( ' ').R().FFP( 1));

      OSTREAM << "  " << count_format( nr_of_residues)
              << count_format( nr_of_chains)
              << count_format( nr_of_ss_bridges)
              << count_format( nr_of_intra_chain_ss_bridges)
              << count_format( nr_of_ss_bridges - nr_of_intra_chain_ss_bridges)
              << "   TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)\n";

      // hydrogenbond summary
      OSTREAM << "  " << count_format( nr_of_hydrogen_bonds)
              << per100_format( nr_of_hydrogen_bonds * 100.0 / nr_of_residues)
              << "   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES\n";

      OSTREAM << "  " << count_format( m_NrOfHydrogenBondsInParallelBridges)
              << per100_format( m_NrOfHydrogenBondsInParallelBridges * 100.0 / nr_of_residues)
              << "   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES\n";

      OSTREAM << "  " << count_format( m_NrOfHydrogenBondsInAntiparallelBridges)
              << per100_format( m_NrOfHydrogenBondsInAntiparallelBridges * 100.0 / nr_of_residues)
              << "   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES\n";

      for( int k( -5); k <= 5; ++k)
      {
        OSTREAM << "  " << count_format( nr_of_hydrogen_bonds_per_distance[ k + 5])
                << per100_format( nr_of_hydrogen_bonds_per_distance[ k + 5] * 100.0 / nr_of_residues)
                << "   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I" << k << "), SAME NUMBER PER 100 RESIDUES5\n";
      }

      // histograms...
      for( size_t i( 1); i <= s_HistogramSize; ++i)
      {
        OSTREAM << count_format( i);
      }
      OSTREAM << "    *** HISTOGRAMS OF ***\n";

      for( size_t i( 0); i < s_HistogramSize; ++i)
      {
        OSTREAM << count_format( nr_res_per_alpha_helix_histogram[ i]);
      }
      OSTREAM << "    RESIDUES PER ALPHA HELIX\n";

      for( size_t i( 0); i < s_HistogramSize; ++i)
      {
        OSTREAM << count_format( m_ParallelBridgesPerLadderHistogram[ i]);
      }
      OSTREAM << "    PARALLEL BRIDGES PER LADDER\n";

      for( size_t i( 0); i < s_HistogramSize; ++i)
      {
        OSTREAM << count_format( m_AntiparallelBridgesPerLadderHistogram[ i]);
      }
      OSTREAM << "    ANTIPARALLEL BRIDGES PER LADDER\n";

      for( size_t i( 0); i < s_HistogramSize; ++i)
      {
        OSTREAM << count_format( m_LaddersPerSheetHistogram[ i]);
      }
      OSTREAM << "    LADDERS PER SHEET\n";

      // per residue information

      OSTREAM << "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA\n";

      // iterate through chains
      for
      (
        assemble::ProteinModel::const_iterator chain_itr( m_Model->GetChains().Begin()), chain_itr_end( m_Model->GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        const util::SiPtrVector< const AABase> res( ( *chain_itr)->GetAminoAcids());
        for
        (
          util::SiPtrVector< const AABase>::const_iterator aa_itr( res.Begin()), aa_itr_end( res.End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          WriteResidue( *aa_itr, OSTREAM);
        }
      }

      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DSSP::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DSSP::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    //! @brief write information for a single residue
    //! @param SP_AMINO_ACID SiPtr to amino acid residue
    //! @param OSTREAM stream to write to
    //! @return the stream written to
    std::ostream &DSSP::WriteResidue( const util::SiPtr< const AABase> &SP_AMINO_ACID, std::ostream &OSTREAM) const
    {
      const Atom &ca( SP_AMINO_ACID->GetAtom( GetAtomTypes().CA));
      const char code( SP_AMINO_ACID->GetType()->GetOneLetterCode());
      char ss( ' ');
      switch( m_SecondaryStructure.find( SP_AMINO_ACID)->second)
      {
        case e_Alphahelix: ss = 'H'; break;
        case e_Betabridge: ss = 'B'; break;
        case e_Strand    : ss = 'E'; break;
        case e_Helix_3   : ss = 'G'; break;
        case e_Helix_5   : ss = 'I'; break;
        case e_Turn      : ss = 'T'; break;
        case e_Bend      : ss = 'S'; break;
        case e_Loop      : ss = ' '; break;
      }

      char helix[ 3];
      for( size_t stride = 3; stride <= 5; ++stride)
      {
        switch( m_HelixFlag[ stride - 3].find( SP_AMINO_ACID)->second)
        {
          case e_HelixNone:        helix[ stride - 3] = ' '; break;
          case e_HelixStart:       helix[ stride - 3] = '>'; break;
          case e_HelixEnd:         helix[ stride - 3] = '<'; break;
          case e_HelixStartAndEnd: helix[ stride - 3] = 'X'; break;
          case e_HelixMiddle:      helix[ stride - 3] = '0' + stride; break;
        }
      }

      const char bend( m_Bend.Contains( SP_AMINO_ACID) ? 'S' : ' ');

      std::pair< double, char> alpha_chirality( 360, ' ');
      if( m_Alpha.find( SP_AMINO_ACID) != m_Alpha.end())
      {
        alpha_chirality = m_Alpha.find( SP_AMINO_ACID)->second;
      }
//      tr1::tie(alpha,chirality) = residue.Alpha();

      size_t bp[ 2] = {};
      char bridgelabel[ 2] = { ' ', ' '};
      for( size_t i( 0); i < 2; ++i)
      {
        const BridgePartner &p( m_BridgePartner.Find( SP_AMINO_ACID)->second( i));
        if( p.m_Partner.IsDefined())
        {
          bp[ i] = p.m_Partner->GetSeqID();
          bp[ i] %= 10000; // won't fit otherwise...
          bridgelabel[ i] = 'A' + p.m_Ladder % 26;
          if( p.m_Parallel)
          {
            bridgelabel[ i] = tolower( bridgelabel[ i]);
          }
        }
      }

      char sheet = ' ';
      if( m_Sheet.Find( SP_AMINO_ACID)->second != 0)
      {
        sheet = 'A' + ( m_Sheet.Find( SP_AMINO_ACID)->second - 1) % 26;
      }

      std::string NHO[ 2], ONH[ 2];
      const storage::VectorND< 2, storage::Pair< util::SiPtr< const AABase>, double> > &acceptors( m_HydrogenBondAcceptor.Find( SP_AMINO_ACID)->second);
      const storage::VectorND< 2, storage::Pair< util::SiPtr< const AABase>, double> > &donors( m_HydrogenBondDonor.Find( SP_AMINO_ACID)->second);
      for( size_t i( 0); i < 2; ++i)
      {
        NHO[ i] = ONH[ i] = "     0, 0.0";

        if( acceptors( i).First().IsDefined())
        {
          const int d( acceptors( i).First()->GetSeqID() - SP_AMINO_ACID->GetSeqID());
          NHO[ i] = util::Format().W( 6).Fill( ' ')( d) + ',' + util::Format().W( 4).FFP( 1)( acceptors( i).Second());
        }

        if( donors( i).First().IsDefined())
        {
          const int d( donors( i).First()->GetSeqID() - SP_AMINO_ACID->GetSeqID());
          ONH[ i] = util::Format().W( 6).Fill( ' ')( d) + ',' + util::Format().W( 4).FFP( 1)( donors( i).Second());
        }
      }

      double kappa( 360);
      if( m_Kappa.find( SP_AMINO_ACID) != m_Kappa.end())
      {
        kappa = m_Kappa.find( SP_AMINO_ACID)->second;
      }

      double phi( 360);
      double psi( 360);
      if( m_Phi.find( SP_AMINO_ACID) != m_Phi.end())
      {
        phi = m_Phi.find( SP_AMINO_ACID)->second;
      }
      if( m_Psi.find( SP_AMINO_ACID) != m_Psi.end())
      {
        psi = m_Psi.find( SP_AMINO_ACID)->second;
      }
      double tco( 0);
      if( m_TCO.find( SP_AMINO_ACID) != m_TCO.end())
      {
        tco = m_TCO.find( SP_AMINO_ACID)->second;
      }

      const util::Format seq_nr_format( util::Format().W( 5).Fill( ' ').R());
      const util::Format bridge_partner_format( util::Format().W( 4).Fill( ' ').R());
      const util::Format coord_format( util::Format().W( 7).Fill( ' ').R().FFP( 1));
      const util::Format angle_format( util::Format().W( 6).Fill( ' ').R().FFP( 1));
      const util::Format tco_format( util::Format().W( 8).Fill( ' ').R().FFP( 3));
      const util::Format acc_format( util::Format().W( 5).Fill( ' ').R());

      OSTREAM << seq_nr_format( SP_AMINO_ACID->GetSeqID())
              << seq_nr_format( SP_AMINO_ACID->GetPdbID())
              << SP_AMINO_ACID->GetPdbICode()
              << SP_AMINO_ACID->GetChainID() << ' '
              << code << "  "
              << ss << ' '
              << helix[ 0] << helix[ 1] << helix[ 2]
              << bend << alpha_chirality.second
              << bridgelabel[ 0] << bridgelabel[ 1]
              << bridge_partner_format( bp[ 0]) << bridge_partner_format( bp[ 1])
              << sheet
              << acc_format( 0)
              << NHO[ 0]
              << ONH[ 0]
              << NHO[ 1]
              << ONH[ 1]
              << tco_format( tco)
              << angle_format( kappa)
              << angle_format( alpha_chirality.first)
              << angle_format( phi)
              << angle_format( psi)
              << coord_format( ca.GetCoordinates().X())
              << coord_format( ca.GetCoordinates().Y())
              << coord_format( ca.GetCoordinates().Z())
              << '\n';

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate all pairwise hydrogen bonding terms
    void DSSP::CalculateHydrogenBondEnergies() const
    {
      const double minimal_ca_distance( 9.0);
      BCL_MessageDbg( "calculate hbond energies");

      const util::SiPtrVector< const AABase> amino_acids( m_Model->GetAminoAcids());

      for( util::SiPtrVector< const AABase>::const_iterator itr_a( amino_acids.Begin()), itr_end( amino_acids.End()); itr_a != itr_end; ++itr_a)
      {
        for( size_t stride( 3); stride <= 5; ++stride)
        {
          m_HelixFlag[ stride - 3][ *itr_a] = e_HelixNone;
        }
        m_BridgePartner[ *itr_a];
        m_SecondaryStructure[ *itr_a] = e_Loop;
        m_Sheet[ *itr_a] = 0;

        for( util::SiPtrVector< const AABase>::const_iterator itr_b( itr_a + 1); itr_b != itr_end; ++itr_b)
        {
          if( Distance( ( *itr_a)->GetAtom( GetAtomTypes().CA), ( *itr_b)->GetAtom( GetAtomTypes().CA)) >= minimal_ca_distance)
          {
            continue;
          }
          CalculateHydrogenBondEnergy( **itr_a, **itr_b);
          if( ( *itr_b)->GetSeqID() != ( *itr_a)->GetSeqID() + 1)
          {
            CalculateHydrogenBondEnergy( **itr_b, **itr_a);
          }
        }

        util::SiPtrVector< const AABase>::const_iterator itr_prev( itr_a);
        util::SiPtrVector< const AABase>::const_iterator itr_current( itr_a + 1);
        if( itr_current == itr_end)
        {
          continue;
        }
        const double phi( ( *itr_current)->CalculatePhi( ( *itr_prev)->GetAtom( GetAtomTypes().C)));
        m_Phi[ *itr_current] = math::Angle::Degree( phi);
        const double psi( ( *itr_prev)->CalculatePsi( ( *itr_current)->GetAtom( GetAtomTypes().N)));
        m_Psi[ *itr_prev] = math::Angle::Degree( psi);

        const double tco
        (
          linal::ProjAngleCosinus
          (
            ( *itr_current)->GetAtom( GetAtomTypes().O).GetCoordinates(),
            ( *itr_current)->GetAtom( GetAtomTypes().C).GetCoordinates(),
            ( *itr_prev)->GetAtom( GetAtomTypes().O).GetCoordinates(),
            ( *itr_prev)->GetAtom( GetAtomTypes().C).GetCoordinates()
          )
        );
        m_TCO[ *itr_current] = tco;

        util::SiPtrVector< const AABase>::const_iterator itr_next( itr_a + 2);
        if( itr_next == itr_end)
        {
          continue;
        }

        util::SiPtrVector< const AABase>::const_iterator itr_next_next( itr_a + 3);
        if( itr_next_next == itr_end)
        {
          continue;
        }

        CalculateAlpha( **itr_prev, **itr_current, **itr_next, **itr_next_next);
      }
    }

    //! @brief using all hydrogen bonds, determine all beta sheets
    void DSSP::CalculateBetaSheets() const
    {
      BCL_MessageDbg( "Calculate beta sheets");
      const util::SiPtrVector< const AABase> amino_acids( m_Model->GetAminoAcids());

      // Calculate Bridges
      std::list< Bridge> bridges;

      if( amino_acids.GetSize() > 4)
      {
        for
        (
          util::SiPtrVector< const AABase>::const_iterator
            itr_i( amino_acids.Begin() + 1), itr_i_end( amino_acids.End() - 4);
          itr_i != itr_i_end;
          ++itr_i
        )
        {
          util::SiPtrVector< const AABase>::const_iterator itr_i_prev( itr_i - 1), itr_i_next( itr_i + 1);

          for
          (
            util::SiPtrVector< const AABase>::const_iterator itr_j( itr_i + 3), itr_j_end( amino_acids.End() - 1);
            itr_j != itr_j_end;
            ++itr_j
          )
          {
            util::SiPtrVector< const AABase>::const_iterator itr_j_prev( itr_j - 1), itr_j_next( itr_j + 1);

            const BridgeType type( TestBridge( **itr_i_prev, **itr_i, **itr_i_next, **itr_j_prev, **itr_j, **itr_j_next));
            if( type == e_BTNoBridge)
            {
              continue;
            }

            bool found( false);
            for
            (
              std::list< Bridge>::iterator bridge_itr( bridges.begin()), bridge_itr_end( bridges.end());
              bridge_itr != bridge_itr_end;
              ++bridge_itr
            )
            {
              Bridge &bridge( *bridge_itr);

              if( type != bridge.m_Type || ( *itr_i)->GetSeqID() != bridge.m_I.back()->GetSeqID() + 1)
              {
                continue;
              }

              if( type == e_BTParallel && bridge.m_J.back()->GetSeqID() + 1 == ( *itr_j)->GetSeqID())
              {
                bridge.m_I.push_back( *itr_i);
                bridge.m_J.push_back( *itr_j);
                found = true;
                break;
              }

              if( type == e_BTAntiParallel && bridge.m_J.front()->GetSeqID() - 1 == ( *itr_j)->GetSeqID())
              {
                bridge.m_I.push_back( *itr_i);
                bridge.m_J.push_front( *itr_j);
                found = true;
                break;
              }
            }

            if( !found)
            {
              Bridge bridge;

              bridge.m_Type = type;
              bridge.m_I.push_back( *itr_i);
              bridge.m_ChainI = ( *itr_i)->GetChainID();
              bridge.m_J.push_back( *itr_j);
              bridge.m_ChainJ = ( *itr_j)->GetChainID();

              bridges.push_back( bridge);
            }
          }
        }
      }

      // extend ladders
      bridges.sort();

      for
      (
        std::list< Bridge>::iterator itr_i( bridges.begin()), itr_end( bridges.end());
        itr_i != itr_end;
        ++itr_i
      )
      {
        std::list< Bridge>::iterator itr_j( itr_i);
        ++itr_j;
        for( ; itr_j != itr_end; ++itr_j)
        {
          const int ibi( itr_i->m_I.front()->GetSeqID());
          const int iei( itr_i->m_I.back()->GetSeqID());
          const int jbi( itr_i->m_J.front()->GetSeqID());
          const int jei( itr_i->m_J.back()->GetSeqID());
          const int ibj( itr_j->m_I.front()->GetSeqID());
          const int iej( itr_j->m_I.back()->GetSeqID());
          const int jbj( itr_j->m_J.front()->GetSeqID());
          const int jej( itr_j->m_J.back()->GetSeqID());

          if
          (
            itr_i->m_Type != itr_j->m_Type ||
//            MResidue::NoChainBreak( inResidues[ min( ibi, ibj)], inResidues[ max( iei, iej)]) == false ||
//            MResidue::NoChainBreak( inResidues[ min( jbi, jbj)], inResidues[ max( jei, jej)]) == false ||
            ibj - iei >= 6 ||
            ( iei >= ibj && ibi <= iej)
          )
          {
            continue;
          }

          bool bulge;
          if( itr_i->m_Type == e_BTParallel)
          {
            bulge = ( ( jbj - jei < 6 && ibj - iei < 3) || ( jbj - jei < 3));
          }
          else
          {
            bulge = ( ( jbi - jej < 6 && ibj - iei < 3) || ( jbi - jej < 3));
          }

          if( bulge)
          {
            itr_i->m_I.insert( itr_i->m_I.end(), itr_j->m_I.begin(), itr_j->m_I.end());
            if( itr_i->m_Type == e_BTParallel)
            {
              itr_i->m_J.insert( itr_i->m_J.end(), itr_j->m_J.begin(), itr_j->m_J.end());
            }
            else
            {
              itr_i->m_J.insert( itr_i->m_J.begin(), itr_j->m_J.begin(), itr_j->m_J.end());
            }
            std::list< Bridge>::iterator itr_erase( itr_j);
            --itr_j;

            bridges.erase( itr_erase);
          }
        }
      }

      // Sheet
      std::list< Bridge *> ladders;
      for( std::list< Bridge>::iterator itr( bridges.begin()), itr_end( bridges.end()); itr != itr_end; ++itr)
      {
        Bridge &bridge( *itr);
        ladders.push_back( &bridge);

        size_t n( bridge.m_I.size());
        if( n > s_HistogramSize)
        {
          n = s_HistogramSize;
        }

        if( bridge.m_Type == e_BTParallel)
        {
          m_ParallelBridgesPerLadderHistogram[ n - 1] += 1;
        }
        else
        {
          m_AntiparallelBridgesPerLadderHistogram[ n - 1] += 1;
        }
      }

      size_t sheet( 1);
      size_t ladder( 0);

      while( !ladders.empty())
      {
        std::list< Bridge *> sheets;
        sheets.splice( sheets.begin(), ladders, ladders.begin());

        bool done( false);
        while( !done)
        {
          done = true;
          for
          (
            std::list< Bridge *>::iterator itr_s( sheets.begin()), itr_s_end( sheets.end());
            itr_s != itr_s_end;
            ++itr_s
          )
          {
            for
            (
              std::list< Bridge *>::iterator itr_l( ladders.begin()), itr_l_end( ladders.end());
              itr_l != itr_l_end;
              ++itr_l
            )
            {
              if( Linked( **itr_s, **itr_l))
              {
                sheets.splice( sheets.end(), ladders, itr_l);
                done = false;
                break;
              }
            }
            if( !done)
            {
              break;
            }
          }
        }

        // create a set of sheets, to store on the bridges
        std::set< Bridge *> sheetset( sheets.begin(), sheets.end());

        for
        (
          std::list< Bridge *>::iterator itr_s( sheets.begin()), itr_s_end( sheets.end());
          itr_s != itr_s_end;
          ++itr_s, ++ladder
        )
        {
          Bridge &bridge( **itr_s);
          bridge.m_Ladder = ladder;
          bridge.m_Sheet = sheet;
          bridge.m_Link = sheetset;
        }

        size_t nr_of_ladders_per_sheet( sheets.size());
        if( nr_of_ladders_per_sheet > s_HistogramSize)
        {
          nr_of_ladders_per_sheet = s_HistogramSize;
        }
        if( nr_of_ladders_per_sheet == 1 && ( *sheets.begin())->m_I.size() > 1)
        {
          m_LaddersPerSheetHistogram[ 0] += 1;
        }
        else if( nr_of_ladders_per_sheet > 1)
        {
          m_LaddersPerSheetHistogram[ nr_of_ladders_per_sheet - 1] += 1;
        }

        ++sheet;
      }

      for( std::list< Bridge>::iterator itr( bridges.begin()), itr_end( bridges.end()); itr != itr_end; ++itr)
      {
        Bridge &bridge( *itr);
        // find out if any of the i and j set members already have
        // a bridge assigned, if so, we're assigning bridge 2

        size_t beta_i( 0);
        size_t beta_j( 0);

        for( std::deque< util::SiPtr< const AABase> >::const_iterator aa_itr( bridge.m_I.begin()), aa_itr_end( bridge.m_I.end()); aa_itr != aa_itr_end; ++aa_itr)
        {
          const BridgePartnerContainerType::const_iterator partner_itr( m_BridgePartner.Find( *aa_itr));
          if( partner_itr != m_BridgePartner.End() && partner_itr->second.First().m_Partner.IsDefined())
          {
            beta_i = 1;
            break;
          }
        }

        for( std::deque< util::SiPtr< const AABase> >::const_iterator aa_itr( bridge.m_J.begin()), aa_itr_end( bridge.m_J.end()); aa_itr != aa_itr_end; ++aa_itr)
        {
          const BridgePartnerContainerType::const_iterator partner_itr( m_BridgePartner.Find( *aa_itr));
          if( partner_itr != m_BridgePartner.End() && partner_itr->second.First().m_Partner.IsDefined())
          {
            beta_j = 1;
            break;
          }
        }

        SecondaryStructureType ss( e_Betabridge);
        if( bridge.m_I.size() > 1)
        {
          ss = e_Strand;
        }

        if( bridge.m_Type == e_BTParallel)
        {
          m_NrOfHydrogenBondsInParallelBridges += bridge.m_I.back()->GetSeqID() - bridge.m_I.front()->GetSeqID() + 2;

          for
          (
            std::deque< util::SiPtr< const AABase> >::const_iterator
              itr_i( bridge.m_I.begin()), itr_j( bridge.m_J.begin()), itr_i_end( bridge.m_I.end());
            itr_i != itr_i_end;
            ++itr_i, ++itr_j
          )
          {
            m_BridgePartner[ *itr_i]( beta_i) = BridgePartner( *itr_j, bridge.m_Ladder, true);
          }

          for
          (
            std::deque< util::SiPtr< const AABase> >::const_iterator
              itr_i( bridge.m_J.begin()), itr_j( bridge.m_I.begin()), itr_i_end( bridge.m_J.end());
            itr_i != itr_i_end;
            ++itr_i, ++itr_j
          )
          {
            m_BridgePartner[ *itr_i]( beta_j) = BridgePartner( *itr_j, bridge.m_Ladder, true);
          }
        }
        else
        {
          m_NrOfHydrogenBondsInAntiparallelBridges += bridge.m_I.back()->GetSeqID() - bridge.m_I.front()->GetSeqID() + 2;

          std::deque< util::SiPtr< const AABase> >::const_reverse_iterator itr_j( bridge.m_J.rbegin());
          for
          (
            std::deque< util::SiPtr< const AABase> >::const_iterator
              itr_i( bridge.m_I.begin()), itr_i_end( bridge.m_I.end());
            itr_i != itr_i_end;
            ++itr_i, ++itr_j
          )
          {
            m_BridgePartner[ *itr_i]( beta_i) = BridgePartner( *itr_j, bridge.m_Ladder, false);
          }

          itr_j = bridge.m_I.rbegin();
          for
          (
            std::deque< util::SiPtr< const AABase> >::const_iterator
              itr_i( bridge.m_J.begin()), itr_i_end( bridge.m_J.end());
            itr_i != itr_i_end;
            ++itr_i, ++itr_j
          )
          {
            m_BridgePartner[ *itr_i]( beta_j) = BridgePartner( *itr_j, bridge.m_Ladder, false);
          }
        }

        for
        (
          std::deque< util::SiPtr< const AABase> >::const_iterator
            aa_itr( bridge.m_I.begin()), aa_itr_end( bridge.m_I.end());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          if( m_SecondaryStructure[ *aa_itr] != e_Strand)
          {
            m_SecondaryStructure[ *aa_itr] = ss;
          }
          m_Sheet[ *aa_itr] = bridge.m_Sheet;
        }

        for
        (
          std::deque< util::SiPtr< const AABase> >::const_iterator
            aa_itr( bridge.m_J.begin()), aa_itr_end( bridge.m_J.end());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          if( m_SecondaryStructure[ *aa_itr] != e_Strand)
          {
            m_SecondaryStructure[ *aa_itr] = ss;
          }
          m_Sheet[ *aa_itr] = bridge.m_Sheet;
        }
      }
    }

    //! @brief using all hydrogen bonds, determine all helix types
    void DSSP::CalculateAlphaHelices() const
    {
      BCL_MessageDbg( "Calculate alpha helices");

      // Helix and Turn
      for
      (
        assemble::ProteinModel::const_iterator chain_itr( m_Model->GetChains().Begin()), chain_itr_end( m_Model->GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        for( size_t stride( 3); stride <= 5; ++stride)
        {
          const util::SiPtrVector< const AABase> res( ( *chain_itr)->GetAminoAcids());
          if( res.GetSize() < stride)
          {
            continue;
          }

          for
          (
            util::SiPtrVector< const AABase>::const_iterator aa_itr( res.Begin()), aa_itr_end( res.End() - stride);
            aa_itr != aa_itr_end;
            ++aa_itr
          )
          {
            if( TestBond( **( aa_itr + stride), **aa_itr)) // and MResidue::NoChainBreak(res[i], res[i + stride]))
            {
              m_HelixFlag[ stride - 3][ *( aa_itr + stride)] = e_HelixEnd;
              for
              (
                util::SiPtrVector< const AABase>::const_iterator aa_j_itr( aa_itr + 1), aa_j_itr_end( aa_itr + stride);
                aa_j_itr != aa_j_itr_end;
                ++aa_j_itr
              )
              {
                if( m_HelixFlag[ stride - 3][ *aa_j_itr] == e_HelixNone)
                {
                  m_HelixFlag[ stride - 3][ *aa_j_itr] = e_HelixMiddle;
                }
              }

              if( m_HelixFlag[ stride - 3][ *aa_itr] == e_HelixEnd)
              {
                m_HelixFlag[ stride - 3][ *aa_itr] = e_HelixStartAndEnd;
              }
              else
              {
                m_HelixFlag[ stride - 3][ *aa_itr] = e_HelixStart;
              }
            }
          }
        }
      }

      // iterate through chains
      for
      (
        assemble::ProteinModel::const_iterator
          chain_itr( m_Model->GetChains().Begin()), chain_itr_end( m_Model->GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        const util::SiPtrVector< const AABase> res( ( *chain_itr)->GetAminoAcids());
        if( res.GetSize() < 5)
        {
          continue;
        }

        for
        (
          util::SiPtrVector< const AABase>::const_iterator aa_itr( res.Begin()), aa_itr_end( res.End() - 4);
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          const AABase &prev_prev( **aa_itr);
          const AABase &center( **( aa_itr + 2));
          const AABase &next_next( **( aa_itr + 4));

          const double kappa( math::Angle::Degree( Kappa( prev_prev, center, next_next)));
          m_Kappa[ util::ToSiPtr( center)] = kappa;
          if( kappa > 70)
          {
            m_Bend.Insert( util::ToSiPtr( center));
          }
        }
      }

      // iterate through chains
      for
      (
        assemble::ProteinModel::const_iterator
          chain_itr( m_Model->GetChains().Begin()), chain_itr_end( m_Model->GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        const util::SiPtrVector< const AABase> res( ( *chain_itr)->GetAminoAcids());
        if( res.GetSize() < 6)
        {
          continue;
        }

        for
        (
          util::SiPtrVector< const AABase>::const_iterator aa_itr( res.Begin()), aa_itr_end( res.End());
          aa_itr + 5 != aa_itr_end;
          ++aa_itr
        )
        {
          if( IsHelixStart( **aa_itr, 4) && IsHelixStart( **( aa_itr + 1), 4))
          {
            for
            (
              util::SiPtrVector< const AABase>::const_iterator aa_j_itr( aa_itr + 1), aa_j_itr_end( aa_j_itr + 4);
              aa_j_itr != aa_j_itr_end;
              ++aa_j_itr
            )
            {
              m_SecondaryStructure[ *aa_j_itr] = e_Alphahelix;
            }
          }
        }
      }

      // iterate through chains
      for
      (
        assemble::ProteinModel::const_iterator
          chain_itr( m_Model->GetChains().Begin()), chain_itr_end( m_Model->GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        const util::SiPtrVector< const AABase> res( ( *chain_itr)->GetAminoAcids());
        if( res.GetSize() < 4)
        {
          continue;
        }

        for
        (
          util::SiPtrVector< const AABase>::const_iterator aa_itr( res.Begin()), aa_itr_end( res.End());
          aa_itr + 4 != aa_itr_end;
          ++aa_itr
        )
        {
          if( IsHelixStart( **aa_itr, 3) && IsHelixStart( **( aa_itr + 1), 3))
          {
            bool empty( true);
            for
            (
              util::SiPtrVector< const AABase>::const_iterator aa_j_itr( aa_itr + 1), aa_j_itr_end( aa_j_itr + 3);
              empty && aa_j_itr != aa_j_itr_end;
              ++aa_j_itr
            )
            {
              empty = m_SecondaryStructure[ *aa_j_itr] == e_Loop || m_SecondaryStructure[ *aa_j_itr] == e_Helix_3;
            }
            if( empty)
            {
              for
              (
                util::SiPtrVector< const AABase>::const_iterator aa_j_itr( aa_itr + 1), aa_j_itr_end( aa_j_itr + 3);
                empty && aa_j_itr != aa_j_itr_end;
                ++aa_j_itr
              )
              {
                m_SecondaryStructure[ *aa_j_itr] = e_Helix_3;
              }
            }
          }
        }
      }

      // iterate through chains
      for
      (
        assemble::ProteinModel::const_iterator chain_itr( m_Model->GetChains().Begin()), chain_itr_end( m_Model->GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        const util::SiPtrVector< const AABase> res( ( *chain_itr)->GetAminoAcids());
        if( res.GetSize() < 6)
        {
          continue;
        }
        for
        (
          util::SiPtrVector< const AABase>::const_iterator aa_itr( res.Begin()), aa_itr_end( res.End());
          aa_itr + 7 != aa_itr_end;
          ++aa_itr
        )
        {
          if( IsHelixStart( **aa_itr, 5) && IsHelixStart( **( aa_itr + 1), 5))
          {
            bool empty( true);
            for
            (
              util::SiPtrVector< const AABase>::const_iterator aa_j_itr( aa_itr + 1), aa_j_itr_end( aa_j_itr + 5);
              empty && aa_j_itr != aa_j_itr_end;
              ++aa_j_itr
            )
            {
              empty = m_SecondaryStructure[ *aa_j_itr] == e_Loop || m_SecondaryStructure[ *aa_j_itr] == e_Helix_5;
            }
            if( empty)
            {
              for
              (
                util::SiPtrVector< const AABase>::const_iterator aa_j_itr( aa_itr + 1), aa_j_itr_end( aa_j_itr + 5);
                empty && aa_j_itr != aa_j_itr_end;
                ++aa_j_itr
              )
              {
                m_SecondaryStructure[ *aa_j_itr] = e_Helix_5;
              }
            }
          }
        }
      }

      // iterate through chains
      for
      (
        assemble::ProteinModel::const_iterator chain_itr( m_Model->GetChains().Begin()), chain_itr_end( m_Model->GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        const util::SiPtrVector< const AABase> res( ( *chain_itr)->GetAminoAcids());
        if( res.GetSize() < 3)
        {
          continue;
        }
        for
        (
          util::SiPtrVector< const AABase>::const_iterator aa_itr( res.Begin() + 1), aa_itr_end( res.End() - 1);
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          if( m_SecondaryStructure[ *aa_itr] != e_Loop)
          {
            continue;
          }
          bool is_turn( false);
          for( size_t stride( 3); stride <= 5 && !is_turn; ++stride)
          {
            for( size_t k( 1); k < stride && !is_turn; ++k)
            {
              is_turn = ( aa_itr >= res.Begin() + k) && IsHelixStart( **( aa_itr - k), stride);
            }
          }

          if( is_turn)
          {
            m_SecondaryStructure[ *aa_itr] = e_Turn;
          }
          else if( m_Bend.Contains( *aa_itr))
          {
            m_SecondaryStructure[ *aa_itr] = e_Bend;
          }
        }
      }
    }

    //! @brief test and calculate hydrogen bond energy between donor and acceptor, store in member
    //! @param AA_DONOR (hydrogen, nitrogen)
    //! @param AA_ACCEPTOR (oxygen, carbon)
    //! @return the energy calculated - if relevant interaction was found
    double DSSP::CalculateHydrogenBondEnergy( const AABase &AA_DONOR, const AABase &AA_ACCEPTOR) const
    {
      const double min_hydrogen_bond_energy( -9.9);
      const double coupling_constant( -27.888); //  = -332 * 0.42 * 0.2
      const double min_distance( 0.5);

      double energy( 0.0);

      const Atom &dh
      (
        AA_DONOR.GetAtom( GetAtomTypes().H).GetCoordinates().IsDefined() ?
            AA_DONOR.GetAtom( GetAtomTypes().H) : AA_DONOR.GetAtom( GetAtomTypes().N)
      );

      if( AA_DONOR.GetType() != GetAATypes().PRO)
      {
        const Atom &dn( AA_DONOR.GetAtom( GetAtomTypes().N));
        const Atom &ao( AA_ACCEPTOR.GetAtom( GetAtomTypes().O));
        const Atom &ac( AA_ACCEPTOR.GetAtom( GetAtomTypes().C));

        const double dist_ho( Distance( dh, ao));
        const double dist_hc( Distance( dh, ac));
        const double dist_nc( Distance( dn, ac));
        const double dist_no( Distance( dn, ao));

        if( dist_ho < min_distance || dist_hc < min_distance || dist_nc < min_distance || dist_no < min_distance)
        {
          energy = min_hydrogen_bond_energy;
        }
        else
        {
          energy = coupling_constant / dist_ho - coupling_constant / dist_hc + coupling_constant / dist_nc - coupling_constant / dist_no;
        }

        // DSSP compatibility mode:
        energy = round( energy * 1000) / 1000;

        if( energy < min_hydrogen_bond_energy)
        {
          energy = min_hydrogen_bond_energy;
        }
      }

      // update donor
      {
        HydrogenBondContainerType::iterator itr( m_HydrogenBondAcceptor.Find( util::ToSiPtr( AA_DONOR)));
        if( itr == m_HydrogenBondAcceptor.End())
        {
          itr = m_HydrogenBondAcceptor.Insert
              (
                std::make_pair
                (
                  util::ToSiPtr( AA_DONOR),
                  storage::VectorND< 2, storage::Pair< util::SiPtr< const AABase>, double> >
                  (
                    storage::Pair< util::SiPtr< const AABase>, double>( util::SiPtr< const AABase>(), 0.0),
                    storage::Pair< util::SiPtr< const AABase>, double>( util::SiPtr< const AABase>(), 0.0)
                  )
                )
              ).first;
        }

        if( energy < itr->second.First().Second())
        {
          itr->second.Second() = itr->second.First();
          itr->second.First() = storage::Pair< util::SiPtr< const AABase>, double>( util::ToSiPtr( AA_ACCEPTOR), energy);
        }
        else if( energy < itr->second.Second().Second())
        {
          itr->second.Second() = storage::Pair< util::SiPtr< const AABase>, double>( util::ToSiPtr( AA_ACCEPTOR), energy);
        }
      }

      // update acceptor
      {
        HydrogenBondContainerType::iterator itr( m_HydrogenBondDonor.Find( util::ToSiPtr( AA_ACCEPTOR)));
        if( itr == m_HydrogenBondDonor.End())
        {
          itr = m_HydrogenBondDonor.Insert
              (
                std::make_pair
                (
                  util::ToSiPtr( AA_ACCEPTOR),
                  storage::VectorND< 2, storage::Pair< util::SiPtr< const AABase>, double> >
                  (
                    storage::Pair< util::SiPtr< const AABase>, double>( util::SiPtr< const AABase>(), 0.0),
                    storage::Pair< util::SiPtr< const AABase>, double>( util::SiPtr< const AABase>(), 0.0)
                  )
                )
              ).first;
        }

        if( energy < itr->second.First().Second())
        {
          itr->second.Second() = itr->second.First();
          itr->second.First() = storage::Pair< util::SiPtr< const AABase>, double>( util::ToSiPtr( AA_DONOR), energy);
        }
        else if( energy < itr->second.Second().Second())
        {
          itr->second.Second() = storage::Pair< util::SiPtr< const AABase>, double>( util::ToSiPtr( AA_DONOR), energy);
        }
      }

      // end
      return energy;
    }

    //! @brief calculate alpha angle, store in member
    //! @param PREV previous amino acid
    //! @param CENTER center amino acid
    //! @param NEXT next amino acid
    //! @param NEXT_NEXT next amino acid after NEXT
    std::pair< double, char> DSSP::CalculateAlpha
    (
      const AABase &PREV, const AABase &CENTER, const AABase &NEXT, const AABase &NEXT_NEXT
    ) const
    {
      const Atom &ca_prev(      PREV.GetAtom(      GetAtomTypes().CA));
      const Atom &ca_current(   CENTER.GetAtom(    GetAtomTypes().CA));
      const Atom &ca_next(      NEXT.GetAtom(      GetAtomTypes().CA));
      const Atom &ca_next_next( NEXT_NEXT.GetAtom( GetAtomTypes().CA));

      const double alpha( math::Angle::Degree( Dihedral( ca_prev, ca_current, ca_next, ca_next_next)));
      char chirality( alpha < 0 ? '-' : '+');

      return m_Alpha[ util::ToSiPtr( CENTER)] = std::make_pair( alpha, chirality);
    }

    //! @brief test if three amino acids in one stretch - I is interacting via a beta bridge with a second stretch - J
    //! @param I_PREV previous amino acid in i
    //! @param I_CURRENT center amino acid in i
    //! @param I_NEXT next amino acid in i
    //! @param J_PREV previous amino acid in j
    //! @param J_CURRENT center amino acid in j
    //! @param J_NEXT next amino acid in j
    //! @return parallel or anti parallel or no bridge was found
    DSSP::BridgeType DSSP::TestBridge
    (
      const AABase &I_PREV, const AABase &I_CURR, const AABase &I_NEXT,
      const AABase &J_PREV, const AABase &J_CURR, const AABase &J_NEXT
    ) const
    {
      // I. a   d II. a   d   parallel
      //      \         /
      //    b   e     b   e
      //      /         \                      ..
      //    c   f     c   f
      //
      // III. a <- f  IV. a   f   antiparallel
      //
      //    b  e      b<->e
      //
      //    c->d      c   d

      BridgeType type( e_BTNoBridge);

//      if( NoChainBreak( A, C) && NoChainBreak( D, F))
      if( I_PREV.GetChainID() == I_NEXT.GetChainID() && J_PREV.GetChainID() == J_NEXT.GetChainID())
      {
        if( ( TestBond( I_NEXT, J_CURR) && TestBond( J_CURR, I_PREV)) || ( TestBond( J_NEXT, I_CURR) && TestBond( I_CURR, J_PREV)))
        {
          type = e_BTParallel;
        }
        else if( ( TestBond( I_NEXT, J_PREV) && TestBond( J_NEXT, I_PREV)) || ( TestBond( J_CURR, I_CURR) && TestBond( I_CURR, J_CURR)))
        {
          type = e_BTAntiParallel;
        }
      }

      return type;
    }

    //! @brief test if donor and acceptor are connected through a hydrogen bond
    //! @param AA_DONOR donor amino acid
    //! @param AA_ACCEPTOR acceptor amino acid
    //! @return true if acceptor is forming a hydrogen bond to donor
    bool DSSP::TestBond( const AABase &AA_DONOR, const AABase &AA_ACCEPTOR) const
    {
      // find the donor entry
      const HydrogenBondContainerType::const_iterator itr( m_HydrogenBondAcceptor.Find( util::ToSiPtr( AA_DONOR)));
      if( itr == m_HydrogenBondAcceptor.End())
      {
        return false;
      }

      return
        ( itr->second.First().First() == util::ToSiPtr( AA_ACCEPTOR) && itr->second.First().Second() < m_MaxHydrogenBondEnergy) ||
        ( itr->second.Second().First() == util::ToSiPtr( AA_ACCEPTOR) && itr->second.Second().Second() < m_MaxHydrogenBondEnergy);
    }

    //! @brief test if any of the residues in bridge a is identical to any of the residues in bridge b in case they
    //!        need to be joined
    //! @param BRIDGE_A bridge a potentially containing any AA of bridge b
    //! @param BRIDGE_B bridge b potentially containing any AA of bridge a
    //! @return true, if bridge and b are overlapping
    bool DSSP::Linked( const Bridge &BRIDGE_A, const Bridge &BRIDGE_B)
    {
      return
        find_first_of( BRIDGE_A.m_I.begin(), BRIDGE_A.m_I.end(), BRIDGE_B.m_I.begin(), BRIDGE_B.m_I.end()) != BRIDGE_A.m_I.end() ||
        find_first_of( BRIDGE_A.m_I.begin(), BRIDGE_A.m_I.end(), BRIDGE_B.m_J.begin(), BRIDGE_B.m_J.end()) != BRIDGE_A.m_I.end() ||
        find_first_of( BRIDGE_A.m_J.begin(), BRIDGE_A.m_J.end(), BRIDGE_B.m_I.begin(), BRIDGE_B.m_I.end()) != BRIDGE_A.m_J.end() ||
        find_first_of( BRIDGE_A.m_J.begin(), BRIDGE_A.m_J.end(), BRIDGE_B.m_J.begin(), BRIDGE_B.m_J.end()) != BRIDGE_A.m_J.end();
    }

    //! @brief calculate kappa angle over 5 residues
    //! @param PREV_PREV center - 2 amino acid
    //! @param CENTER center amino acid
    //! @param NEXT_NEXT center + 2 amino acid
    //! @return the kappa angle
    double DSSP::Kappa( const AABase &PREV_PREV, const AABase &CENTER, const AABase &NEXT_NEXT)
    {
      double result( 2 * math::g_Pi);

      const double ckap
      (
        linal::ProjAngleCosinus
        (
          NEXT_NEXT.GetCA().GetCoordinates(),
          CENTER.GetCA().GetCoordinates(),
          CENTER.GetCA().GetCoordinates(),
          PREV_PREV.GetCA().GetCoordinates()
        )
      );

      const double skap( sqrt( 1.0 - ckap * ckap));

      result = atan2( skap, ckap);

      return result;
    }

    //! @brief test if given amino acid is a start residue for given helix stride
    //! @param AMINO_ACID the amino acid in question
    //! @param STRIDE the stride of the helix
    //! @return true if amino acid was flagged start compatible for the stride
    bool DSSP::IsHelixStart( const AABase &AMINO_ACID, const size_t STRIDE) const
    {
      HelixFlagContainerType::const_iterator itr( m_HelixFlag[ STRIDE - 3].find( util::ToSiPtr( AMINO_ACID)));
      if( itr == m_HelixFlag[ STRIDE - 3].end())
      {
        return false;
      }
      return itr->second == e_HelixStart || itr->second == e_HelixStartAndEnd;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_environment_type_data.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_environment_types.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> EnvironmentTypeData::s_Instance( GetObjectInstances().AddInstance( new EnvironmentTypeData()));

    //! @brief construct undefined environment type
    EnvironmentTypeData::EnvironmentTypeData() :
      m_TwoLetterCode( ""),
      m_ReducedType( ""),
      m_ReducedIndex( util::GetUndefined< size_t>()),
      m_IsGap( false),
      m_DefaultThickness( util::GetUndefined< double>())
    {
    }

    //! @brief construct EnvironmentTypeData
    //! @param TWO_LETTER_CODE two letter code for this environment type
    //! @param REDUCED_TYPE_NAME two letter code of reduced EnvironmentType
    //! @param REDUCED_INDEX index for the reduced EnvironmentType
    //! @param IS_GAP boolean to indicate whether this type corresponds to a gap region
    //! @param DEFAULT_THICKNESS efault thickness for this environment type
    EnvironmentTypeData::EnvironmentTypeData
    (
      const std::string &TWO_LETTER_CODE,
      const std::string &REDUCED_TYPE_NAME,
      const size_t REDUCED_INDEX,
      const bool IS_GAP,
      const double DEFAULT_THICKNESS
    ) :
      m_TwoLetterCode( TWO_LETTER_CODE),
      m_ReducedType( REDUCED_TYPE_NAME),
      m_ReducedIndex( REDUCED_INDEX),
      m_IsGap( IS_GAP),
      m_DefaultThickness( DEFAULT_THICKNESS)
    {
    }

    //! @brief virtual copy constructor
    EnvironmentTypeData *EnvironmentTypeData::Clone() const
    {
      return new EnvironmentTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &EnvironmentTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get two letter code
    //! @return two letter code
    const std::string &EnvironmentTypeData::GetTwoLetterCode() const
    {
      return m_TwoLetterCode;
    }

    //! @brief get index for reduced type
    //! @return index for reduced type
    size_t EnvironmentTypeData::GetReducedIndex() const
    {
      return m_ReducedIndex;
    }

    //! @brief get the reduced environment type
    //! this maps environment types to one of three states: Core, Transition or Solution
    //! @return reduced EnvironmentType
    const EnvironmentType &EnvironmentTypeData::GetReducedType() const
    {
      return GetEnvironmentTypes().EnvironmentTypeFromTwoLetterCode( m_ReducedType);
    }

    //! @brief get the reduced EnvironmentType's name
    //! @return name of reduced EnvironmentType for this type
    const std::string &EnvironmentTypeData::GetReducedTypeString() const
    {
      return m_ReducedType;
    }

    //! @brief get default thickness
    double EnvironmentTypeData::GetDefaultThickness() const
    {
      return m_DefaultThickness;
    }

    //! @brief get whether this is a gap region
    //! @return whether this is a gap region
    bool EnvironmentTypeData::IsGap() const
    {
      return m_IsGap;
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &EnvironmentTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_TwoLetterCode, ISTREAM);
      io::Serialize::Read( m_ReducedType, ISTREAM);
      io::Serialize::Read( m_ReducedIndex, ISTREAM);
      io::Serialize::Read( m_IsGap, ISTREAM);
      io::Serialize::Read( m_DefaultThickness, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &EnvironmentTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_TwoLetterCode, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ReducedType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ReducedIndex, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IsGap, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DefaultThickness, OSTREAM, INDENT) << '\n';

      // return
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_environment_types.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief construct all EnvironmentTypes
    EnvironmentTypes::EnvironmentTypes() :
      util::Enumerate< EnvironmentTypeData, EnvironmentTypes>( false),
      e_MembraneCore(          AddEnum( "MEMBRANE_CORE",                    EnvironmentTypeData( "MC", "MC", 0, false, double(     10.0)))),
      e_GapCoreTransition(     AddEnum( "MEMBRANE_GAP_CORE_TRANSITION",     EnvironmentTypeData( "G1", "TR", 1, true,  double(      2.5)))),
      e_Transition(            AddEnum( "MEMBRANE_TRANSITION",              EnvironmentTypeData( "TR", "TR", 1, false, double(     10.0)))),
      e_GapTransitionSolution( AddEnum( "MEMBRANE_GAP_TRANSITION_SOLUTION", EnvironmentTypeData( "G2", "SO", 2, true,  double(      2.5)))),
      e_Solution(              AddEnum( "SOLUTION",                         EnvironmentTypeData( "SO", "SO", 2, false, double( 100000.0)))),
      e_MembraneInside(        AddEnum( "MEMBRANE_INSIDE",                  EnvironmentTypeData( "MI", "MC", 0, false, double(     10.0)))),
      e_MembraneOutside(       AddEnum( "MEMBRANE_OUTSIDE",                 EnvironmentTypeData( "ME", "MC", 0, false, double(     10.0)))),
      e_SolutionInside(        AddEnum( "SOLUTION_INSIDE",                  EnvironmentTypeData( "SI", "SO", 2, false, double( 100000.0)))),
      e_SolutionOutside(       AddEnum( "SOLUTION_OUTSIDE",                 EnvironmentTypeData( "SE", "SO", 2, false, double( 100000.0))))
    {
    }

    //! @brief function to deduce EnvironmentType from two letter code
    //! @param TWO_LETTER_CODE two letter code descriptor for environment type of interest
    //! @return EnvironmentType specified by the given TWO_LETTER_CODE
    const EnvironmentType &EnvironmentTypes::EnvironmentTypeFromTwoLetterCode( const std::string &TWO_LETTER_CODE) const
    {
      // iterate over all EnvironmentTypes
      for
      (
        const_iterator type_itr( Begin()), type_itr_end( End());
        type_itr != type_itr_end;
        ++type_itr
      )
      {
        // if ONE_LETTER_CODE matches this EnvironmentTypes's two letter code
        if( ( *type_itr)->GetTwoLetterCode() == TWO_LETTER_CODE)
        {
          // return this EnvironmentType
          return *type_itr;
        }
      }

      // if no match was found return undefined EnvironmentType
      return GetEnvironmentTypes().e_Undefined;
    }

    //! @brief returns number of reduced types
    //! @return number of reduced types
    size_t EnvironmentTypes::GetNumberReducedTypes() const
    {
      // end
      return GetReducedTypes().GetSize();
    }

    //! @brief returns a vector that contains the 3 reduces types, e_MembraneCore, e_Transition and e_Solution
    //! @return a vector that contains the 3 reduces types, e_MembraneCore, e_Transition and e_Solution
    const storage::Vector< EnvironmentType> &EnvironmentTypes::GetReducedTypes() const
    {
      // initialize a static vector to store reduced types
      static const storage::Vector< EnvironmentType> s_reduced_types_vector
      (
        storage::Vector< EnvironmentType>::Create( e_MembraneCore, e_Transition, e_Solution)
      );

      // end
      return s_reduced_types_vector;
    }

    //! @brief construct on access function for all EnvironmentTypes
    //! @return reference to only instance of EnvironmentTypes enum
    const EnvironmentTypes &GetEnvironmentTypes()
    {
      return EnvironmentTypes::GetEnums();
    }

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< biol::EnvironmentTypeData, biol::EnvironmentTypes>;

  } // namespace util
} // namespace bcl
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
#include "biol/bcl_biol_exposure_prediction.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_running_average.h"
#include "model/bcl_model_interface_retrieve_from_file.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ExposurePrediction::s_Instance
    (
      GetObjectInstances().AddInstance( new ExposurePrediction())
    );

    //! @brief return model path as a string
    //! @param ExposureType one of the exposure types
    //! @return model path as a string
    const std::string &ExposurePrediction::GetModelPath( const ExposureType &EXPOSURE_TYPE)
    {
      static std::string s_paths[] =
      {
          "exposure/contact_number/protomeric",
          "exposure/contact_number/oligomeric",
          "exposure/rsa/protomeric",
          "exposure/rsa/oligomeric",
          GetStaticClassName< ExposurePrediction::ExposureType>()
      };
      return s_paths[ size_t( EXPOSURE_TYPE)];
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &ExposurePrediction::GetFileExtension( const ExposureType &EXPOSURE_TYPE)
    {
      static const std::string s_extensions[] =
      {
          ".proto_cn",
          ".oligo_cn",
          ".proto_rsa",
          ".oligo_rsa",
          GetStaticClassName< ExposurePrediction::ExposureType>()
      };
      return s_extensions[ size_t( EXPOSURE_TYPE)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ExposurePrediction::ExposurePrediction()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ExposurePrediction
    ExposurePrediction *ExposurePrediction::Clone() const
    {
      return new ExposurePrediction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ExposurePrediction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief iterates over the sequences and calculates the jufo predictions for every residue in the sequence
    //! @param SEQUENCE sequence of interest
    //! @param EXPOSURE_TYPE one of the exposure types
    void ExposurePrediction::Calculate( AASequence &SEQUENCE, const ExposureType &EXPOSURE_TYPE)
    {
      // create a dummy model for the prediction
      assemble::ProteinModel model
      (
        util::ShPtr< assemble::Chain>
        (
          new assemble::Chain
          (
            util::ShPtr< AASequence>( new AASequence( SEQUENCE))
          )
        )
      );
      ExposurePrediction::Calculate( model, EXPOSURE_TYPE);
    }

    //! @brief iterates over the sequences in ProteinModel and calculates the jufo predictions for every residue in the sequence
    //! @param PROTEIN_MODEL ProteinModel for which JUFO will be calculated
    //! @param EXPOSURE_TYPE one of the exposure types
    void ExposurePrediction::Calculate( assemble::ProteinModel &PROTEIN_MODEL, const ExposureType &EXPOSURE_TYPE)
    {
      // create a protein-model-with-cache
      assemble::ProteinModelWithCache pmwc( PROTEIN_MODEL, false);
      Calculate( pmwc, EXPOSURE_TYPE);
    }

    //! @brief iterates over the sequences in ProteinModel and calculates the jufo predictions for every residue in the sequence
    //! @param PROTEIN_MODEL ProteinModel for which JUFO will be calculated
    //! @param EXPOSURE_TYPE one of the exposure types
    void ExposurePrediction::Calculate( assemble::ProteinModelWithCache &PROTEIN_MODEL, const ExposureType &EXPOSURE_TYPE)
    {
      // make sure the blast profile exists for this sequence
      BCL_Assert
      (
        PROTEIN_MODEL.GetIterator()->GetBlastProfilePtr().IsDefined(),
        "Blast Profile is not available"
      );

      // path to the model
      const std::string path_to_model( GetModelPath( EXPOSURE_TYPE));

      // create the descriptor to generate the dataset
      util::Implementation< descriptor::Base< AABase, float> > aa_descriptor
      (
        "PredictionMean(storage=File(directory=" + model::Model::AddModelPath( path_to_model) + ", prefix=model))"
      );

      // set the object up
      aa_descriptor->SetObject( PROTEIN_MODEL);

      // set the dimension (1 because we operate on elements of the sequence)
      aa_descriptor->SetDimension( 1);

      // create a descriptor iterator
      descriptor::Iterator< AABase> itr( PROTEIN_MODEL.GetIterator());

      // create a non-const iterator over the same amino acids
      iterate::Generic< AABase> itr_non_const( PROTEIN_MODEL.GetIteratorNonConst());

      // iterate over the amino acids
      for( ; itr.NotAtEnd(); ++itr, ++itr_non_const)
      {
        // set the prediction
        itr_non_const->SetExposurePrediction( aa_descriptor->operator ()( itr)( 0));
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads the exposure predictions for a model given a path and prefix
    //! @param PROTEIN_MODEL protein model to contain predictions
    //! @param PREFIX prefix of the sequence ( usually pdb id)
    //! @param PATH path where exposure prediction files can be found
    void ExposurePrediction::ReadPredictions
    (
      assemble::ProteinModel &PROTEIN_MODEL,
      const std::string &PREFIX,
      const std::string &PATH,
      const ExposureType &EXPOSURE_TYPE
    )
    {
      //iterate over all chains in the sequence
      for
      (
        util::ShPtrVector< assemble::Chain>::iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // construct the filename
        std::string filename
        (
          PATH + PATH_SEPARATOR + PREFIX + ( *chain_itr)->GetChainID() + GetFileExtension( EXPOSURE_TYPE)
        );

        if( !io::DirectoryEntry( filename).DoesExist())
        {
          filename = PATH + PATH_SEPARATOR + PREFIX + GetFileExtension( EXPOSURE_TYPE);
        }
        // open the file
        io::IFStream read;
        io::File::MustOpenIFStream( read, filename);

        // read predictions for the sequence
        ReadPredictions( read, *( *chain_itr)->GetSequence());
        io::File::CloseClearFStream( read);
      }
    }

    //! @brief reads the exposure predictions from a file
    //! @param ISTREAM stream to read from
    //! @param PROTEIN_MODEL protein model to contain predictions
    void ExposurePrediction::ReadPredictions( std::istream &ISTREAM, assemble::ProteinModel &PROTEIN_MODEL)
    {
      //iterate over all chains in the sequence
      for
      (
        util::ShPtrVector< assemble::Chain>::iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // read predictions for the sequence
        ReadPredictions( ISTREAM, *( *chain_itr)->GetSequence());
      }
    }

    //! @brief reads the exposure predictions from a file
    //! @param ISTREAM stream to read from
    //! @param SEQUENCE sequence to contain predictions
    void ExposurePrediction::ReadPredictions( std::istream &ISTREAM, AASequence &SEQUENCE)
    {
      // iterate over the sequence
      for
      (
        AASequence::iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // read in the seq id and exposure
        int seq_id;
        std::string exposure;
        ISTREAM >> seq_id >> exposure;

        BCL_Assert
        (
          ( *aa_itr)->GetSeqID() == seq_id,
          "Mismatch in seq id for exposure predictions. " + ( *aa_itr)->GetIdentification() + " does not match " +
            util::Format()( seq_id)
        );

        // set the prediction
        ( *aa_itr)->SetExposurePrediction( util::ConvertStringToNumericalValue< double>( exposure));
      }
    }

    //! @brief writes out the exposure predictions to a file
    //! @param OSTREAM stream to write to
    //! @param SEQUENCE sequence containing predictions
    void ExposurePrediction::WritePredictions( std::ostream &OSTREAM, const AASequence &SEQUENCE)
    {
      // iterate over the sequence
      for
      (
        AASequence::const_iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        OSTREAM << ( *aa_itr)->GetSeqID() << '\t' << ( *aa_itr)->GetExposurePrediction() << '\n';
      }
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ExposurePrediction::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ExposurePrediction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_membrane.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_environment_types.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! @brief return command line flag for defining the membrane thickness
    //! @return command line flag for defining the membrane thickness
    util::ShPtr< command::FlagInterface> &Membrane::GetFlagMembrane()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic( "membrane", "flag for using a membrane and defining membrane thicknesses")
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters
        flag->PushBack( GetParameterCoreThickness());
        flag->PushBack( GetParameterTransitionThickness());
        flag->PushBack( GetParameterGapThickness());
      }

      // end
      return s_flag;
    }

    //! @brief return command line parameter for defining the membrane core thickness
    //! @return command line parameter for defining the membrane core thickness
    util::ShPtr< command::ParameterInterface> &Membrane::GetParameterCoreThickness()
    {
      // initialize static instance of the parameter
      static util::ShPtr< command::ParameterInterface> s_core_thickness_param
      (
        new command::Parameter
        (
          "membrane_core",
          "half the thickness of membrane core region in Angstroem",
          command::ParameterCheckRanged< double>( 0.0, 50.0),
          "10.0"
        )
      );

      // end
      return s_core_thickness_param;
    }

    //! @brief return command line parameter for defining the membrane transition thickness
    //! @return command line parameter for defining the membrane transition thickness
    util::ShPtr< command::ParameterInterface> &Membrane::GetParameterTransitionThickness()
    {
      // initialize static instance of the parameter
      static util::ShPtr< command::ParameterInterface> s_transition_thickness_param
      (
        new command::Parameter
        (
          "membrane_trans",
          "thickness of membrane transition region in Angstroem",
          command::ParameterCheckRanged< double>( 0.0, 50.0),
          "10.0"
        )
      );

      // end
      return s_transition_thickness_param;
    }

    //! @brief return command line parameter for defining the membrane gap thickness
    //! @return command line parameter for defining the membrane gap thickness
    util::ShPtr< command::ParameterInterface> &Membrane::GetParameterGapThickness()
    {
      // initialize static instance of the parameter
      static util::ShPtr< command::ParameterInterface> s_gap_thickness_param
      (
        new command::Parameter
        (
          "membrane_gap",
          "thickness of membrane gap between regions in Angstroem",
          command::ParameterCheckRanged< double>( 0.0, 50.0),
          "2.5"
        )
      );

      // end
      return s_gap_thickness_param;
    }

    //! @brief create a membrane object from commandline arguments
    //! @return membrane object created from commandline arguments
    Membrane Membrane::GetCommandLineMembrane()
    {
      // construct and return a membrane object constructed from commandline arguments
      return Membrane
      (
        GetParameterCoreThickness()->GetNumericalValue< double>(),
        GetParameterTransitionThickness()->GetNumericalValue< double>(),
        GetParameterGapThickness()->GetNumericalValue< double>()
      );
    }

    //! @brief static undefined membrane object
    //! @return undefined membrane object
    const Membrane &Membrane::GetUndefinedMembrane()
    {
      static const Membrane s_membrane
      (
        util::GetUndefined< double>(),
        util::GetUndefined< double>(),
        util::GetUndefined< double>()
      );
      return s_membrane;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Membrane::s_Instance
    (
      GetObjectInstances().AddInstance( new Membrane())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Membrane::Membrane() :
      m_Orientation(),
      m_Thicknesses
      (
        FillThicknessVector
        (
          GetEnvironmentTypes().e_MembraneCore->GetDefaultThickness(),
          GetEnvironmentTypes().e_Transition->GetDefaultThickness(),
          GetEnvironmentTypes().e_GapCoreTransition->GetDefaultThickness()
        )
      ),
      m_Limits( FillLimitsVector( m_Thicknesses))
    {
      m_IsDefined = m_Orientation.IsDefined() && IsDefined( m_Thicknesses) && IsDefined( m_Limits);
    }

    //! @brief constructor from membrane normal, all thicknesses and gap thickness
    //! @param THICKNESSES is a vector of membrane thicknesses
    //! @param CENTER center of the membrane
    //! @param NORMAL membrane normal
    Membrane::Membrane
    (
      const storage::Vector< double> &THICKNESSES,
      const linal::Vector3D &NORMAL,
      const linal::Vector3D &CENTER
    ) :
      m_Orientation( OrientationFromNormal( NORMAL, CENTER)),
      m_Thicknesses( THICKNESSES),
      m_Limits( FillLimitsVector( m_Thicknesses))
    {
      m_IsDefined = m_Orientation.IsDefined() && IsDefined( m_Thicknesses) && IsDefined( m_Limits);
      BCL_Assert
      (
        THICKNESSES.GetSize() <= GetEnvironmentTypes().GetEnumCount(),
        "given vector of thicknesses had too many elements"
      );
    }

    //! @brief constructor from membrane normal, specified thicknesses and gap thickness
    //! @param THICKNESS_CORE thickness of membrane core
    //! @param THICKNESS_TRANSITION thickness of membrane transition region
    //! @param THICKNESS_GAP thickness of membrane gap
    //! @param NORMAL Membrane normal
    Membrane::Membrane
    (
      const double THICKNESS_CORE,
      const double THICKNESS_TRANSITION,
      const double THICKNESS_GAP,
      const linal::Vector3D &NORMAL
    ) :
      m_Orientation( OrientationFromNormal( NORMAL)),
      m_Thicknesses( FillThicknessVector( THICKNESS_CORE, THICKNESS_TRANSITION, THICKNESS_GAP)),
      m_Limits( FillLimitsVector( m_Thicknesses))
    {
      m_IsDefined = m_Orientation.IsDefined() && IsDefined( m_Thicknesses) && IsDefined( m_Limits);
    }

    //! @brief virtual copy constructor
    Membrane *Membrane::Clone() const
    {
      return new Membrane( *this);
    }

    //! @brief destructor
    Membrane::~Membrane()
    {
      m_DestructorSignal.Emit( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Membrane::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the membrane normal
    //! @return the membrane normal
    linal::Vector3D Membrane::GetNormal() const
    {
      linal::Vector3D normal( GetAxis( coord::GetAxes().e_Z));
      normal.Normalize();
      return normal;
    }

    //! @brief return the Thickness of the region
    //! @return the Thickness of the region
    double Membrane::GetThickness( const EnvironmentType &ENVIRONMENT) const
    {
      return m_Thicknesses( ENVIRONMENT);
    }

    //! @brief return the Limit of the region
    //! @return the Limit of the region
    double Membrane::GetLimit( const EnvironmentType &ENVIRONMENT) const
    {
      return m_Limits( ENVIRONMENT);
    }

    //! @brief returns true if the membrane is not set to static undefined membrane
    //! @return true if the membrane is not set to static undefined membrane
    bool Membrane::IsDefined() const
    {
      return m_IsDefined;
    }

    //! @brief returns the geometric center of the object
    //! @return the geometric center of the object
    linal::Vector3D Membrane::GetCenter() const
    {
      return m_Orientation.GetOrigin();
    }

    //! @brief return the orientation of the object
    //! @return orientation
    linal::Vector3D Membrane::GetAxis( const coord::Axis &AXIS) const
    {
      return m_Orientation.GetAxis( AXIS);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns EnvironmentType according to the coordinates
    //! @param COORDINATES coordinates of object
    //! @return EnvironmentType according to the coordinates
    EnvironmentType Membrane::DetermineEnvironmentType( const linal::Vector3D &COORDINATES) const
    {
      // return undefined environment for undefined coords
      if( !COORDINATES.IsDefined())
      {
        return GetEnvironmentTypes().e_Undefined;
      }

      // calculate the distance from the membrane plane
      const double distance( DistanceFromPlane( COORDINATES));

      //check in which region the distance is the first time smaller than the current limit
      for
      (
        EnvironmentTypes::const_iterator itr( GetEnvironmentTypes().Begin()), itr_end( GetEnvironmentTypes().End());
        itr != itr_end; ++itr
      )
      {
        if( distance < m_Limits( *itr))
        {
          return *itr;
        }
      }

      //if it does not fit in any region, return undefined
      return util::GetUndefined< EnvironmentType>();
    }

    //! @brief returns EnvironmentType according to the z-coordiante and weight, if it is in a gap reagion
    //! @param COORDINATES coordinates of object
    //! @return pair of EnvironmentType and a weight between 0 and one, 1 if it is closer to the innermose region, 0 if it is close to the outer most region
    storage::Pair< EnvironmentType, double>
    Membrane::DetermineEnvironmentTypeAndWeight( const linal::Vector3D &COORDINATES) const
    {
      //determine the absolute z-coordinate
      const double distance( DistanceFromPlane( COORDINATES));

      storage::Pair< EnvironmentType, double> environment_weight( DetermineEnvironmentType( COORDINATES), 1.0);
      const EnvironmentType env_type( environment_weight.First());

      // switch over environment types and change the weight, if a gap gets hit
      if( env_type == GetEnvironmentTypes().e_GapCoreTransition)
      {
        const double position_in_gap
        (
          ( distance - m_Limits( GetEnvironmentTypes().e_MembraneCore)) /
            m_Thicknesses( GetEnvironmentTypes().e_GapCoreTransition) * math::g_Pi
        );
        environment_weight.Second() = math::WeightBetweenZeroAndPi( position_in_gap);
      }
      else if( env_type == GetEnvironmentTypes().e_GapTransitionSolution)
      {
        const double position_in_gap
        (
          ( distance - m_Limits( GetEnvironmentTypes().e_Transition)) /
            m_Thicknesses( GetEnvironmentTypes().e_GapTransitionSolution) * math::g_Pi);
        environment_weight.Second() = math::WeightBetweenZeroAndPi( position_in_gap);
      }
      // if undefined
      else if( !env_type.IsDefined())
      {
        environment_weight.Second() = util::GetUndefined< double>();
      }

      // return determined environment and weight
      return environment_weight;
    }

    //! @brief  returns solvation energy for given z-coordinate and a vector containing three state solvation energies
    //! @param COORDINATES coordinates of object
    //! @return solvation energy for given z-coordinate and a vector containing three state solvation energies
    double Membrane::CalculateSolvationEnergy( const linal::Vector3D &COORDINATES, const linal::Vector3D &TFE) const
    {
      // determine environment and weight
      const storage::Pair< EnvironmentType, double> env_weight( DetermineEnvironmentTypeAndWeight( COORDINATES));

      // if env type is undefined
      if( !env_weight.First().IsDefined())
      {
        return util::GetUndefined< double>();
      }

      // for actual region
      if( !env_weight.First()->IsGap())
      {
        return TFE( env_weight.First()->GetReducedIndex());
      }

      // for gap
      return           env_weight.Second()  * TFE( env_weight.First()->GetReducedIndex() - 1)
             + ( 1.0 - env_weight.Second()) * TFE( env_weight.First()->GetReducedIndex());
    }

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION Translation to be applied
    void Membrane::Translate( const linal::Vector3D &TRANSLATION)
    {
      m_Orientation( TRANSLATION);
    }

    //! @brief transform the object by a given TransformationMatrix3D
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void Membrane::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      m_Orientation( TRANSFORMATION_MATRIX_3D);
    }

    //! @brief rotate the object by a given RotationMatrix3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
    void Membrane::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      m_Orientation( ROTATION_MATRIX_3D);
    }

    //! @brief fills the Thicknesses with given values
    //! @param THICKNESS_CORE thickness of membrane core region
    //! @param THICKNESS_TRANSITION thickness of membrane transition region
    //! @param THICKNESS_GAP thickness of membrane gap region
    //! @return Thickness vector
    storage::Vector< double> Membrane::FillThicknessVector
    (
      const double THICKNESS_CORE,
      const double THICKNESS_TRANSITION,
      const double THICKNESS_GAP
    )
    {
      storage::Vector< double> thicknesses( GetEnvironmentTypes().GetEnumCount());
      thicknesses( GetEnvironmentTypes().e_MembraneCore)          = THICKNESS_CORE;
      thicknesses( GetEnvironmentTypes().e_GapCoreTransition)     = THICKNESS_GAP;
      thicknesses( GetEnvironmentTypes().e_Transition)            = THICKNESS_TRANSITION;
      thicknesses( GetEnvironmentTypes().e_GapTransitionSolution) = THICKNESS_GAP;
      thicknesses( GetEnvironmentTypes().e_Solution)              = GetEnvironmentTypes().e_Solution->GetDefaultThickness();
      thicknesses( GetEnvironmentTypes().e_SolutionInside)        = 0.0;
      thicknesses( GetEnvironmentTypes().e_SolutionOutside)       = 0.0;

      return thicknesses;
    }

    //! @brief Calculates the Limits for the membrane regions
    //! @param THICKNESSES vector of thickness for membrane regions
    //! @return limits for the membrane regions
    storage::Vector< double> Membrane::FillLimitsVector( const storage::Vector< double> &THICKNESSES)
    {
      // vector of limits
      storage::Vector< double> limits( GetEnvironmentTypes().GetEnumCount());

      // set the value
      double sum_thickness( 0.0);

      // for all other limits the limit is the previous limit + the current thickness
      for
      (
        EnvironmentTypes::const_iterator itr( GetEnvironmentTypes().Begin()),
          itr_end( GetEnvironmentTypes().End());
        itr != itr_end; ++itr
      )
      {
        if( itr->GetIndex() < THICKNESSES.GetSize())
        {
          sum_thickness += THICKNESSES( *itr);
        }
        limits( *itr) = sum_thickness;
      }

      return limits;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Membrane::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Orientation, ISTREAM);
      io::Serialize::Read( m_Thicknesses, ISTREAM);
      io::Serialize::Read( m_Limits, ISTREAM);

      //end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &Membrane::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      //write member
      io::Serialize::Write( m_Orientation, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Thicknesses, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Limits, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    math::TransformationMatrix3D TransformationMatrixFromPDBTMXML_TMATRIX( std::istream &ISTREAM)
    {
      // create linebuffer
      std::string line_buffer;

      linal::Matrix< double> transformation( 4, 4, double( 0.0));

      std::getline( ISTREAM, line_buffer);
      storage::Vector< std::string> xyz( util::SplitString( line_buffer, "="));
      BCL_Assert( util::TrimString( xyz( 0)) == "<ROWX X", "improper format");
      transformation( 0, 0) = util::ConvertStringToNumericalValue< double>( xyz( 1).substr( 1, xyz( 1).length() - 4));
      transformation( 1, 0) = util::ConvertStringToNumericalValue< double>( xyz( 2).substr( 1, xyz( 2).length() - 4));
      transformation( 2, 0) = util::ConvertStringToNumericalValue< double>( xyz( 3).substr( 1, xyz( 3).length() - 4));
      transformation( 3, 0) = util::ConvertStringToNumericalValue< double>( xyz( 4).substr( 1, xyz( 4).length() - 4));
      std::getline( ISTREAM, line_buffer);
      xyz = util::SplitString( line_buffer, "=");
      BCL_Assert( util::TrimString( xyz( 0)) == "<ROWY X", "improper format");
      transformation( 0, 1) = util::ConvertStringToNumericalValue< double>( xyz( 1).substr( 1, xyz( 1).length() - 4));
      transformation( 1, 1) = util::ConvertStringToNumericalValue< double>( xyz( 2).substr( 1, xyz( 2).length() - 4));
      transformation( 2, 1) = util::ConvertStringToNumericalValue< double>( xyz( 3).substr( 1, xyz( 3).length() - 4));
      transformation( 3, 1) = util::ConvertStringToNumericalValue< double>( xyz( 4).substr( 1, xyz( 4).length() - 4));
      std::getline( ISTREAM, line_buffer);
      xyz = util::SplitString( line_buffer, "=");
      BCL_Assert( util::TrimString( xyz( 0)) == "<ROWZ X", "improper format");
      transformation( 0, 2) = util::ConvertStringToNumericalValue< double>( xyz( 1).substr( 1, xyz( 1).length() - 4));
      transformation( 1, 2) = util::ConvertStringToNumericalValue< double>( xyz( 2).substr( 1, xyz( 2).length() - 4));
      transformation( 2, 2) = util::ConvertStringToNumericalValue< double>( xyz( 3).substr( 1, xyz( 3).length() - 4));
      transformation( 3, 2) = util::ConvertStringToNumericalValue< double>( xyz( 4).substr( 1, xyz( 4).length() - 4));

      // check for valid end
      std::getline( ISTREAM, line_buffer);
      BCL_Assert( util::TrimString( line_buffer) == "</TMATRIX>", "improper format");

      // end
      return math::TransformationMatrix3D( transformation);
    }

    //! @brief create membrane object and transformation matrix from given pdbtm xml file
    //! only the core thickness can be retrieved from the xml file, so that transition region and gap are passed
    //! @param ISTREAM input stream
    //! @param THICKNESS_TRANSITION thickness of the membrane transition region
    //! @param THICKNESS_GAP thickness of the gaps between the regions
    //! @return pair of membrane and transformation matrix
    storage::Pair< Membrane, math::TransformationMatrix3D> Membrane::MembraneAndTransformationFromPDBTMXML
    (
      std::istream &ISTREAM,
      const double THICKNESS_TRANSITION,
      const double THICKNESS_GAP
    )
    {
      // create membrane and transformation matrix, initialize with undefined values
      storage::Pair< Membrane, math::TransformationMatrix3D> membrane_matrix
      (
        Membrane( storage::Vector< double>( GetEnvironmentTypes().GetEnumCount(), util::GetUndefined< double>()), linal::Vector3D()),
        math::TransformationMatrix3D( util::UndefinedObject())
      );

      // create linebuffer
      std::string line_buffer;
      while( std::getline( ISTREAM, line_buffer))
      {
        if( util::TrimString( line_buffer) == "<MEMBRANE>")
        {
          break;
        }
      }

      // if no membrane tag was found, then presumably the protein was soluble, return an undefined membrane
      if( !ISTREAM.good())
      {
        return membrane_matrix;
      }

      // get membrane normale
      std::getline( ISTREAM, line_buffer);
      storage::Vector< std::string> xyz( util::SplitString( line_buffer, "="));
      linal::Vector3D normale;
      BCL_Assert( util::TrimString( xyz( 0)) == "<NORMAL X", "improper format, expected \"<NORMAL X\" but found " + util::TrimString( xyz( 0)));
      normale.X() = util::ConvertStringToNumericalValue< double>( xyz( 1).substr( 1, xyz( 1).length() - 4));
      normale.Y() = util::ConvertStringToNumericalValue< double>( xyz( 2).substr( 1, xyz( 2).length() - 4));
      normale.Z() = util::ConvertStringToNumericalValue< double>( xyz( 3).substr( 1, xyz( 3).length() - 4));

      // get and generate matrix
      std::getline( ISTREAM, line_buffer);
      BCL_Assert( util::TrimString( line_buffer) == "<TMATRIX>", "improper format");

      // prepare
      membrane_matrix.First() = Membrane( normale.Norm(), THICKNESS_TRANSITION, THICKNESS_GAP, normale);
      membrane_matrix.Second() = TransformationMatrixFromPDBTMXML_TMATRIX( ISTREAM);

      // check for valid end
      std::getline( ISTREAM, line_buffer);
      BCL_Assert( util::TrimString( line_buffer) == "</MEMBRANE>", "improper format");

      // some debug information
      BCL_MessageDbg( "transformation matrix found in pdbtm xml file:\n" + util::Format()( membrane_matrix.Second()));

      // end
      return membrane_matrix;
    }

    //! @brief returns pairs of chain id and transformation matrix that has to be applied to the chain to get one monomer for
    //! the biological relevant unit
    //! @param ISTREAM input stream of pdbtm xml file
    //! @param CHAIN_IDS chain ids in the original pdb
    //! @return vector of chain id the transformation is applied to, the new chain id and the matrix that needs to be applied order to get the relevant BIOMOLECULE
    storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> >
    Membrane::BioTransformationMatricesFromPDBTMXML( std::istream &ISTREAM, const std::string &CHAIN_IDS)
    {
      std::string delete_chain_ids;
      storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > bio_transformations;

      // create linebuffer and progress to BIOMATRIX
      std::string line_buffer;
      while( std::getline( ISTREAM, line_buffer) && util::TrimString( line_buffer) != "<BIOMATRIX>");

      while( util::TrimString( line_buffer) != "</BIOMATRIX>")
      {
        // error in getting next line
        if( !std::getline( ISTREAM, line_buffer).good())
        {
          break;
        }
        // split line
        storage::Vector< std::string> split_identifier( util::SplitString( line_buffer, "="));

        // check if chain delete was given
        if( util::TrimString( split_identifier( 0)) == "<DELETE CHAINID")
        {
          delete_chain_ids.push_back( split_identifier( 1)[ 1]);
          continue;
        }

        // check if their is a new matrix
        if( util::TrimString( split_identifier( 0)) == "<MATRIX ID")
        {
          // gather all chain ids that this needs to be applied to with its new chain id
          storage::Vector< storage::Pair< char, char> > chain_ids;

          std::getline( ISTREAM, line_buffer);
          // split line
          storage::Vector< std::string> split_apply( util::SplitString( line_buffer, "="));
          while( util::TrimString( split_apply( 0)) == "<APPLY_TO_CHAIN CHAINID")
          {
            const storage::Pair< char, char> chain_id_pair( split_apply( 1)[ 1], split_apply( 2)[ 1]);

            // if the chain id is present in the model
            if( CHAIN_IDS.find( chain_id_pair.First()) != std::string::npos)
            {
              // add the pair
              chain_ids.PushBack( chain_id_pair);
            }
            std::getline( ISTREAM, line_buffer);
            split_apply = util::SplitString( line_buffer, "=");
          }

          // get the matrix for those chains
          BCL_Assert
          (
            util::TrimString( split_apply( 0)) == "<TMATRIX>",
            "improper format. expected <TMATRIX> found: " + split_apply( 0)
          );
          const math::TransformationMatrix3D current_transformation( TransformationMatrixFromPDBTMXML_TMATRIX( ISTREAM));
          std::getline( ISTREAM, line_buffer);
          BCL_Assert
          (
            util::TrimString( line_buffer) == "</MATRIX>",
            "improper format. expected </MATRIX> found: " + line_buffer
          );

          for( storage::Vector< storage::Pair< char, char> >::const_iterator itr( chain_ids.Begin()), itr_end( chain_ids.End()); itr != itr_end; ++itr)
          {
            bio_transformations.PushBack( storage::Triplet< char, char, math::TransformationMatrix3D>( itr->First(), itr->Second(), current_transformation));
          }
        }
      }

      // no bio transformations were found and no chain to delete
      if( bio_transformations.IsEmpty() && delete_chain_ids.empty())
      {
        return bio_transformations;
      }

      // add identity to all given chains, since they are not given in the pdbtmxml file
      for( std::string::const_iterator itr( CHAIN_IDS.begin()), itr_end( CHAIN_IDS.end()); itr != itr_end; ++itr)
      {
        // identify matrix; old chain id is new chain id
        bio_transformations.PushBack( storage::Triplet< char, char, math::TransformationMatrix3D>( *itr, *itr, math::TransformationMatrix3D()));
      }

      // delete undesired chains
      for( std::string::const_iterator char_itr( delete_chain_ids.begin()), char_itr_end( delete_chain_ids.end()); char_itr != char_itr_end; ++char_itr)
      {
        for
        (
          storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> >::iterator
            itr( bio_transformations.Begin()), itr_end( bio_transformations.End()); itr != itr_end; ++itr
        )
        {
          if( itr->Second() == *char_itr)
          {
            bio_transformations.Remove( itr);
            break;
          }
        }
      }

      // end
      return bio_transformations;
    }

    //! @brief evaluates whether all entries in the vector are defined
    //! @param VECTOR vector to be evaulated
    //! @return whether all entries in the vector are defined
    bool Membrane::IsDefined( const storage::Vector< double> &VECTOR)
    {
      // iterate through the vector
      for
      (
        storage::Vector< double>::const_iterator itr( VECTOR.Begin()), itr_end( VECTOR.End());
        itr != itr_end; ++itr
      )
      {
        // if the value is not defined
        if( !util::IsDefined( *itr))
        {
          return false;
        }
      }

      // if this point is reached all values are defined
      return true;
    }

    //! @brief calculates the distance of a point from the center of the membrane plane
    //! @param COORDINATES coordinates to measure
    //! @return calculated distance
    double Membrane::DistanceFromPlane( const linal::Vector3D &COORDINATES) const
    {
      return math::Absolute
      (
        linal::ScalarProduct( GetNormal(), COORDINATES - GetCenter())
      );
    }

    //! @brief gets the membrane orientation from the given normal
    //! @param NORMAL membrane normal
    //! @param CENTER membrane center
    //! @return membrane orientation
    math::TransformationMatrix3D Membrane::OrientationFromNormal
    (
      const linal::Vector3D &NORMAL,
      const linal::Vector3D &CENTER
    )
    {
      // initialize transformation
      math::TransformationMatrix3D transform;

      // get the axis and angle
      linal::Vector3D normal( NORMAL);
      normal.Normalize();
      const linal::Vector3D axis( linal::CrossProduct( normal, coord::GetAxes().e_Z));
      const double angle( ProjAngle( normal, coord::GetAxes().e_Z));

      // only change the orientation if the axis and angle are not zero
      if( axis.Norm() != 0.0 && angle != 0.0)
      {
        transform( math::RotationMatrix3D( axis, angle));
      }

      // apply translation
      transform( CENTER);

      // end
      return transform;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_mutation.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "biol/bcl_biol_atom.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_string_functions.h"
#include "util/bcl_util_string_numeric_conversion.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new Mutation
    Mutation *Mutation::Clone() const
    {
      return new Mutation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Mutation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief write the mutation as a string in standard mutation format
    std::string Mutation::ToString() const
    {
      return std::string( size_t( 1), m_NativeType->GetOneLetterCode())
             + util::Format()( m_ResidueNumber)
             + std::string( size_t( 1), m_MutantType->GetOneLetterCode());
    }

    //! @brief write the mutation as a string in standard mutation format
    Mutation Mutation::FromString( const std::string &STR)
    {
      BCL_Assert( STR.length() >= size_t( 3), "need at least 3 characters for a valid mutation!");
      size_t end( STR.find_first_not_of( "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"));
      if( end == std::string::npos)
      {
        end = STR.size();
      }
      return
        Mutation
        (
          util::ConvertStringToNumericalValue< int>( STR.substr( 1, end - 2)),
          GetAATypes().AATypeFromOneLetterCode( STR[0]),
          GetAATypes().AATypeFromOneLetterCode( STR[end -1]),
          end + size_t( 1) < STR.size() ? STR.substr( end + size_t( 1)) : ""
        );
    }

    //! @param WITH_DATA whether to include any data members, else, only include initialization members
    util::ObjectDataLabel Mutation::GetLabel( const bool &WITH_DATA) const
    {
      return util::ObjectDataLabel( "", ToString());
    }

    //! @brief connect the mutation with a particular residue
    void Mutation::SetAA( const AABase &BASE) const
    {
      for( auto itr( m_AAs.Begin()), itr_end( m_AAs.End()); itr != itr_end; ++itr)
      {
        if( ( *itr)->GetSeqID() == BASE.GetSeqID() && ( *itr)->GetChainID() == BASE.GetChainID())
        {
          *itr = util::ToSiPtr( BASE);
          return;
        }
      }
      m_AAs.PushBack( util::ToSiPtr( BASE));
      if( m_ChainIDs.find( BASE.GetChainID()) == std::string::npos)
      {
        m_ChainIDs += BASE.GetChainID();
      }
      BCL_MessageVrb( "New aa: " + BASE.GetIdentification() + " " + util::Format()( BASE.GetCA().GetCoordinates()));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief set the value of the corresponding member based on the label
    //! @param LABEL label that is used to set the string
    //! @param ERROR_STREAM stream to write errors to
    //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
    bool Mutation::TryRead( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      if( LABEL.GetValue().size() < size_t( 3))
      {
        this->WriteHelp( ERROR_STREAM);
        return false;
      }
      m_NativeType = GetAATypes().AATypeFromOneLetterCode( LABEL.GetValue()[ 0]);
      if( !m_NativeType.IsDefined())
      {
        ERROR_STREAM << LABEL.GetValue()[ 0] << " in mutation " << LABEL.ToString() << " is not a valid amino acid type";
        return false;
      }
      size_t end( LABEL.GetValue().find_first_not_of( "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"));
      if( end == std::string::npos)
      {
        end = LABEL.GetValue().size();
      }
      m_MutantType = GetAATypes().AATypeFromOneLetterCode( LABEL.GetValue()[ end - size_t( 1)]);
      if( !m_MutantType.IsDefined())
      {
        ERROR_STREAM << LABEL.GetValue()[ LABEL.GetValue().size() - size_t( 1)] << " in mutation "
                     << LABEL.ToString() << " is not a valid amino acid type";
        return false;
      }
      if( !util::TryConvertFromString( m_ResidueNumber, LABEL.GetValue().substr( 1, end - 2), ERROR_STREAM))
      {
        ERROR_STREAM << " in mutation " << LABEL.ToString();
        return false;
      }
      m_ChainIDs = end + size_t( 1) < LABEL.GetValue().size() ? LABEL.GetValue().substr( end + size_t( 1)) : "";
      return true;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Mutation::Read( std::istream &ISTREAM)
    {
      BCL_Assert( this->TryRead( util::ObjectDataLabel( ISTREAM), util::GetLogger()), "Could not read mutation");
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &Mutation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM << ToString();
    }

    //! @brief writes the help for the label
    //! @param OSTREAM the stream to which the help is written to
    //! @param INDENT the amount of indent
    //! @return the given stream to which the help was written to
    std::ostream &Mutation::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM << "Mutation in standard format: NativeAATypeResIDMutantAAType like V215M. May optionally indicate chain"
                     << " id(s) afterwards, e.g. V215M-ABCD";
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_protein_charge.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_protein_params.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
  //////////
  // data //
  //////////

    const AATypeData::PropertyType ProteinCharge::s_FirstpKAAProperty   = AATypeData::e_pK_EMBOSS;
    const AATypeData::PropertyType ProteinCharge::s_LastpKAAProperty    = AATypeData::e_pK_ProMoST;

    const double ProteinCharge::s_pK_NTerm[ s_LastpKAAProperty - s_FirstpKAAProperty + 1] =
    {
      8.6, 8.0, 9.6, 8.2, 8.0, 11.2, 8.2, 9.69, 7.7, 0.0, 0.0
    };

    const double ProteinCharge::s_pK_CTerm[ s_LastpKAAProperty - s_FirstpKAAProperty + 1] =
    {
      3.6, 3.1, 2.4, 3.2, 3.1, 4.2, 3.65, 2.34, 3.3, 0.0, 0.0
    };

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinCharge::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinCharge())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinCharge::ProteinCharge() :
      m_AACount(),
      m_NTermAA(),
      m_CTermAA(),
      m_pKProperty( s_LastpKAAProperty)
    {
    }

    //! @brief constructor from amino acid sequence
    //! @param SEQUENCE the amino acid seuqence of the protein of interest
    ProteinCharge::ProteinCharge( const AASequence &SEQUENCE) :
      m_AACount( ProteinParams::CountAAs( SEQUENCE)),
      m_pKProperty( s_LastpKAAProperty)
    {
      if( SEQUENCE.GetSize() < 1)
      {
        return;
      }
      m_NTermAA = SEQUENCE.GetFirstAA()->GetType();
      m_CTermAA = SEQUENCE.GetLastAA()->GetType();
    }

    //! @brief construct from amino acid count
    //! @param AA_COUNT the number of amino acids
    //! @param N_TERM_AA n terminal amino acid
    //! @param C_TERM_AA c terminal amino acid
    ProteinCharge::ProteinCharge
    (
      const storage::Map< AAType, size_t> &AA_COUNT,
      const AAType &N_TERM_AA,
      const AAType &C_TERM_AA
    ) :
      m_AACount( AA_COUNT),
      m_NTermAA( N_TERM_AA),
      m_CTermAA( C_TERM_AA),
      m_pKProperty( s_LastpKAAProperty)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProteinCharge
    ProteinCharge *ProteinCharge::Clone() const
    {
      return new ProteinCharge( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinCharge::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinCharge::GetScheme() const
    {
      return this->GetClassIdentifier();
    }

    //! @brief set the pk aa property to be used
    //! @brief AA_PROPERTY the aa property to be used; needs to be within the allowed pk property range
    void ProteinCharge::SetPKProperty( const AATypeData::PropertyType &PROPERTY)
    {
      if( ( PROPERTY >= s_FirstpKAAProperty && PROPERTY <= s_LastpKAAProperty))
      {
        m_pKProperty = PROPERTY;
      }
      else
      {
        BCL_MessageStd
        (
          "setting the pk property to the default property: " +
          AATypeData::PropertyTypeEnum( s_LastpKAAProperty).GetString()
        );
        m_pKProperty = s_LastpKAAProperty;
      }
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief for a given pH value clauclate the charge
    //! @param PH the pH of the solution
    //! @return the net charge of the protein within this pH
    double ProteinCharge::operator()( const double &PH) const
    {
      double net_charge( 0.0);

      for
      (
        storage::Map< AAType, size_t>::const_iterator itr( m_AACount.Begin()), itr_end( m_AACount.End());
        itr != itr_end;
        ++itr
      )
      {
        net_charge += double( itr->second) * HendersonHasselbachEquation( PH, itr->first->GetAAProperty( m_pKProperty), itr->first->GetAAProperty( AATypeData::e_Charge));
      }

      net_charge = CorrectTerminalCharge( net_charge, PH);

      return net_charge;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinCharge::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AACount, ISTREAM);
      io::Serialize::Read( m_NTermAA, ISTREAM);
      io::Serialize::Read( m_CTermAA, ISTREAM);
      AATypeData::PropertyTypeEnum property;
      io::Serialize::Read( property, ISTREAM);
      SetPKProperty( property);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinCharge::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AACount, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NTermAA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CTermAA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( AATypeData::PropertyTypeEnum( m_pKProperty), OSTREAM, INDENT) << '\n';

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief correct the given charge by the terminal residue charges for that ph
    //! @param PH the desired pH
    //! @param CHARGE ignoring the termini
    //! @return the correct charge by the terminal amino acids, that usually have different pK value
    double ProteinCharge::CorrectTerminalCharge( const double &CHARGE, const double &PH) const
    {
      double corrected_charge( CHARGE);

      // correct termini for ProMoST
      if( m_pKProperty == AATypeData::e_pK_ProMoST)
      {
//        corrected_charge -= HendersonHasselbachEquation( PH, m_NTermAA->GetAAProperty( m_pKProperty), m_NTermAA->GetAAProperty( AATypeData::e_Charge));
//        corrected_charge -= HendersonHasselbachEquation( PH, m_CTermAA->GetAAProperty( m_pKProperty), m_CTermAA->GetAAProperty( AATypeData::e_Charge));
        corrected_charge += HendersonHasselbachEquation( PH, m_NTermAA->GetAAProperty( AATypeData::e_pK_ProMoST_NTerm),  1.0);
        corrected_charge += HendersonHasselbachEquation( PH, m_CTermAA->GetAAProperty( AATypeData::e_pK_ProMoST_CTerm), -1.0);
      }

      // correct termini for Bjellqvist
      else if( m_pKProperty == AATypeData::e_pK_Bjellqvist)
      {
//        corrected_charge -= HendersonHasselbachEquation( PH, m_NTermAA->GetAAProperty( m_pKProperty), m_NTermAA->GetAAProperty( AATypeData::e_Charge));
//        corrected_charge -= HendersonHasselbachEquation( PH, m_CTermAA->GetAAProperty( m_pKProperty), m_CTermAA->GetAAProperty( AATypeData::e_Charge));
        corrected_charge += HendersonHasselbachEquation( PH, m_NTermAA->GetAAProperty( AATypeData::e_pK_Bjellqvist_NTerm),  1.0);
        corrected_charge += HendersonHasselbachEquation( PH, m_CTermAA->GetAAProperty( AATypeData::e_pK_Bjellqvist_CTerm), -1.0);
      }

      // for all other pK vlaue sources, C and N termini are the same for every amino acid
      else
      {
        corrected_charge += HendersonHasselbachEquation( PH, s_pK_CTerm[ m_pKProperty - s_FirstpKAAProperty],  1.0);
        corrected_charge += HendersonHasselbachEquation( PH, s_pK_NTerm[ m_pKProperty - s_FirstpKAAProperty], -1.0);
      }

      return corrected_charge;
    }

    //! @brief Henderson-Hasselbach-Equation
    //! @param PH the ph
    //! @param PK the pK
    //! @param AA_CHARGE the charge of the amino acid (positive or negative)
    //! @return the net charge - if pk is 0 charge returned is 0
    double ProteinCharge::HendersonHasselbachEquation( const double &PH, const double &PK, const double &AA_CHARGE)
    {
      if( PK == 0.0 || AA_CHARGE == 0.0)
      {
        return 0;
      }

      if( AA_CHARGE < 0.0)
      {
        return -1.0 / ( 1.0 + math::Pow( 10.0, PK - PH));
      }

      return 1.0 / ( 1.0 + math::Pow( 10.0, PH - PK));
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_protein_mutation_set.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ProteinMutationSet::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinMutationSet())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinMutationSet::ProteinMutationSet()
    {
    }

    //! @brief construct from a protein model
    //! @param MODEL protein model of interest
    //! @param REQUIRE_COORDINATES whether to exclude sspred analysis methods from AAs that lac defined coordinates
    ProteinMutationSet::ProteinMutationSet
    (
      const assemble::ProteinModel &MODEL,
      const bool &REQUIRE_COORDINATES,
      const storage::Vector< Mutation> &POSSIBLE_MUTATIONS
    ) :
      descriptor::SequenceInterface< Mutation>(),
      m_MutantModel( assemble::ProteinModel::HardCopy( MODEL), REQUIRE_COORDINATES),
      m_PossibleMutations( POSSIBLE_MUTATIONS)
    {
    }

    //! @brief copy constructor
    //! @param ORIGINAL model with cache to copy
    ProteinMutationSet::ProteinMutationSet( const ProteinMutationSet &ORIGINAL) :
      descriptor::SequenceInterface< Mutation>( ORIGINAL),
      m_MutantModel( assemble::ProteinModel::HardCopy( ORIGINAL.m_MutantModel), ORIGINAL.m_MutantModel.GetRequiresCoordinates()),
      m_PossibleMutations( ORIGINAL.m_PossibleMutations)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer a new ProteinMutationSet copied from this model
    ProteinMutationSet *ProteinMutationSet::Clone() const
    {
      return new ProteinMutationSet( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinMutationSet::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the length of the sequence in question
    //! @return the length of the sequence in question
    size_t ProteinMutationSet::GetSize() const
    {
      return m_PossibleMutations.GetSize();
    }

    //! @brief get the iterator for the sequence
    //! @return the iterator for the sequence
    iterate::Generic< const Mutation> ProteinMutationSet::GetIterator() const
    {
      return iterate::Generic< const Mutation>( m_PossibleMutations.Begin(), m_PossibleMutations.End());
    }

    //! @brief get a non-constant iterator for the sequence
    //! @return the non-constant iterator for the sequence
    iterate::Generic< Mutation> ProteinMutationSet::GetIteratorNonConst()
    {
      return iterate::Generic< Mutation>( m_PossibleMutations.Begin(), m_PossibleMutations.End());
    }

    //! @brief Reset the cache
    void ProteinMutationSet::ResetCache() const
    {
      descriptor::SequenceInterface< Mutation>::ResetCache();
      m_MutantModel.ResetCache();
    }

    //! @brief get a particular mutant protein model
    //! @return the mutated protein model
    const assemble::ProteinModelWithMutations &ProteinMutationSet::GetMutant( const Mutation &MUTATION) const
    {
      if( m_MutantModel.OnlyHasMutation( MUTATION))
      {
        return m_MutantModel;
      }
      m_MutantModel.RevertToWildType();
      m_MutantModel.Mutate( MUTATION);
      return m_MutantModel;
    }

    //! @brief get a particular mutant protein model
    //! @return the mutated protein model
    const assemble::ProteinModelWithMutations &ProteinMutationSet::GetNativeType() const
    {
      if( !m_MutantModel.IsWildType())
      {
        m_MutantModel.RevertToWildType();
      }
      return m_MutantModel;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator
    //! @param PROTEIN_MODEL ProteinMutationSet to be copied
    //! @return This model after all members are assigned to values from PROTEIN_MODEL
    ProteinMutationSet &ProteinMutationSet::operator =( const ProteinMutationSet &PROTEIN_MODEL)
    {
      // update members
      m_MutantModel = PROTEIN_MODEL.m_MutantModel;
      m_PossibleMutations = PROTEIN_MODEL.m_PossibleMutations;

      // end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read ProteinMutationSet from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinMutationSet::Read( std::istream &ISTREAM)
    {
      //read data
      io::Serialize::Read( m_MutantModel, ISTREAM);
      io::Serialize::Read( m_PossibleMutations, ISTREAM);

      //return
      return ISTREAM;
    }

    //! @brief write ProteinMutationSet to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &ProteinMutationSet::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write protein model
      io::Serialize::Write( m_MutantModel, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PossibleMutations, OSTREAM, INDENT);
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_protein_params.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_protein_charge.h"
#include "opti/bcl_opti_approximator_root_regula_falsi.h"
#include "opti/bcl_opti_criterion_convergence_argument.h"
#include "storage/bcl_storage_table.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinParams::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinParams())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinParams::ProteinParams()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProteinParams
    ProteinParams *ProteinParams::Clone() const
    {
      return new ProteinParams( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinParams::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate protein parameters
    //! @param SEQUENCE amino acid sequence
    //! @return table containing all information
    storage::Table< double> ProteinParams::operator()( const AASequence &SEQUENCE) const
    {
      static const double s_weight_water
      (
        2 * chemistry::GetElementTypes().e_Hydrogen->GetProperty( chemistry::ElementTypeData::e_Mass) +
        1 * chemistry::GetElementTypes().e_Oxygen->GetProperty( chemistry::ElementTypeData::e_Mass)
      );

      // result
      storage::Table< double> result( storage::TableHeader( storage::Vector< std::string>::Create( "value", "percentage")));

      // count the number of amino acids of each type
      const storage::Map< AAType, size_t> number_aas( CountAAs( SEQUENCE));

      // add the aa count
      AddAACountToResultTable( number_aas, result);

      // molecular weight
      storage::Row< double> &mol_weight_row( result.InsertRow( "molecular_weight[Da]"));
      mol_weight_row( 0) = CalcualteMolecularWeight( number_aas) + s_weight_water;
      mol_weight_row( 1) = 100;

      // extinction coefficient
      CalcualteExtinctionCoefficient( number_aas, result);

      if( SEQUENCE.GetSize() < 1)
      {
        return result;
      }

      // pI
      for
      (
        AATypeData::PropertyTypeEnum pk_property( ProteinCharge::s_FirstpKAAProperty);
        pk_property <= ProteinCharge::s_LastpKAAProperty;
        ++pk_property
      )
      {
        storage::Row< double> &pi_row( result.InsertRow( "pI_" + pk_property.GetString()));
        pi_row( 0) = CalculatePI( number_aas, SEQUENCE.GetFirstAA()->GetType(), SEQUENCE.GetLastAA()->GetType(), pk_property);
        pi_row( 1) = 100;
      }

      // end
      return result;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinParams::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinParams::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the molecular weight from the number of amino acids
    //! @param AA_COUNT a count for each amino acid type
    //! @return molecular weight in Dalton [Da]
    double ProteinParams::CalcualteMolecularWeight( const storage::Map< AAType, size_t> &AA_COUNT)
    {
      double mol_weight( 0.0);

      for
      (
        storage::Map< AAType, size_t>::const_iterator itr( AA_COUNT.Begin()), itr_end( AA_COUNT.End());
        itr != itr_end;
        ++itr
      )
      {
        mol_weight += itr->second * itr->first->GetAAProperty( AATypeData::e_Mass);
      }

      return mol_weight;
    }

    //! @brief calculate the pI value (the pH where the protein is not charged)
    //! @param AA_COUNT a count for each amino acid type
    //! @param N_TERM_TYPE type of n terminal amino acid
    //! @param C_TERM_TYPE type of c terminal amino acid
    //! @param PK_PROPERTY the pk value scale to use
    //! @return the PI value - the pH at which the protein would have no net-charge
    double ProteinParams::CalculatePI
    (
      const storage::Map< AAType, size_t> &AA_COUNT,
      const AAType &N_TERM_AATYPE,
      const AAType &C_TERM_AATYPE,
      const AATypeData::PropertyType &PK_PROPERTY
    )
    {
      // borders of the interval - ph 0 to 14
      const double border_left( 0);
      const double border_right( 14);

      // tolerance for the approximation result
      const double root_tolerance( 0.00001);

      // create termination criteria
      opti::CriterionConvergenceArgument< double, double> criterion( 1, root_tolerance);

      util::ShPtr< ProteinCharge> protein_charge_function( new ProteinCharge( AA_COUNT, N_TERM_AATYPE, C_TERM_AATYPE));
      protein_charge_function->SetPKProperty( PK_PROPERTY);

      // create the approximator
      opti::ApproximatorRootRegulaFalsi< double, double> approximator
      (
        *protein_charge_function, criterion, border_left, border_right
      );

      // approximation
      approximator.Approximate();
      const util::ShPtr< storage::Pair< double, double> > sp_result( approximator.GetTracker().GetBest());

      return sp_result->First();
    }

    //! @brief calculate the extinction coefficient from the number of amino acids in units of  M-1 cm-1, at 280 nm measured in water
    //! @param AA_COUNT a count for each amino acid type
    //! @param RESULT result table
    void ProteinParams::CalcualteExtinctionCoefficient
    (
      const storage::Map< AAType, size_t> &AA_COUNT, storage::Table< double> &RESULT
    )
    {
      const storage::Map< AAType, double> coefficient_map( ExtinctionCoefficientMap());

      double coefficient( 0);
      double coefficient_cystines( 0);
      for
      (
        storage::Map< AAType, double>::const_iterator itr( coefficient_map.Begin()), itr_end( coefficient_map.End());
        itr != itr_end;
        ++itr
      )
      {
        if( AA_COUNT.Has( itr->first))
        {
          if( itr->first == GetAATypes().CYS)
          {
            // only even number of cysteins considered, since only peptide bond does contribute
            coefficient_cystines += ( AA_COUNT.Find( itr->first)->second / 2) * ( 2 * itr->second);
            continue;
          }
          else
          {
            coefficient += AA_COUNT.Find( itr->first)->second * itr->second;
          }
        }
      }

      RESULT.InsertRow( "ExtinctionCoefficient[M-1*cm-1]")( 0) = coefficient;
      RESULT.InsertRow( "ExtinctionCoefficientCystines[M-1*cm-1]")( 0) = coefficient + coefficient_cystines;
    }

    //! @brief count the number of each amino acid in the given sequence
    //! @param SEQUENCE amino acid sequence
    //! @return the number of each amino acid type
    storage::Map< AAType, size_t> ProteinParams::CountAAs( const AASequence &SEQUENCE)
    {
      // count the number of amino acids of each type
      storage::Map< AAType, size_t> number_aas;

      for( AASequence::const_iterator itr( SEQUENCE.Begin()), itr_end( SEQUENCE.End()); itr != itr_end; ++itr)
      {
        ++number_aas[ ( *itr)->GetType()];
      }

      // end
      return number_aas;
    }

    //! @brief extinction coefficient map at 280 nm in water
    //! @return map of amino acid extinction coefficients
    storage::Map< AAType, double> ProteinParams::ExtinctionCoefficientMap()
    {
      storage::Map< AAType, double> map;

      map[ GetAATypes().TRP] = 5500.0;
      map[ GetAATypes().TYR] = 1490.0;
      map[ GetAATypes().CYS] =   62.5; // 125 for each cystine - which are two sulifde-bonded cysteins

      return map;
    }

    //! @brief add the aa count to the result table
    //! @param AA_COUNT map of aatype and their counts
    //! @param RESULT result table
    void ProteinParams::AddAACountToResultTable
    (
      const storage::Map< AAType, size_t> &AA_COUNT,
      storage::Table< double> &RESULT
    )
    {
      storage::Map< AAType, size_t> aa_count( AA_COUNT);
      size_t number_aas( 0);

      // iterate through all natural aa types
      for
      (
        AATypes::const_iterator itr( GetAATypes().Begin()), itr_end( GetAATypes().VAL.GetIterator() + 1);
        itr != itr_end;
        ++itr
      )
      {
        const AAType &current_type( *itr);
        const std::string name( current_type->GetThreeLetterCode() + '_' + current_type->GetOneLetterCode());
        storage::Row< double> &new_row( RESULT.InsertRow( name));
        new_row( 0) = aa_count[ current_type];
        number_aas += new_row( 0);
        aa_count.Erase( current_type);
      }

      // add the count of non-natural aa types
      for
      (
        storage::Map< AAType, size_t>::const_iterator itr( aa_count.Begin()), itr_end( aa_count.End());
        itr != itr_end;
        ++itr
      )
      {
        const AAType &current_type( itr->first);
        const std::string name( current_type->GetThreeLetterCode() + '_' + current_type->GetOneLetterCode());
        storage::Row< double> &new_row( RESULT.InsertRow( name));
        new_row( 0) = itr->second;
        number_aas += itr->second;
      }

      // total
      storage::Row< double> &new_row( RESULT.InsertRow( "number_aas"));
      new_row( 0) = number_aas;

      // percentage
      for
      (
        storage::Table< double>::iterator itr( RESULT.Begin()), itr_end( RESULT.End());
        itr != itr_end && itr->First() != "number_aas";
        ++itr
      )
      {
        itr->Second()( 1) = size_t( ( itr->Second()( 0) * 1000) / number_aas) / 10.0;
      }

      // neg charged residues
      {
        size_t neg_count( 0);
        if( AA_COUNT.Has( GetAATypes().ASP))
        {
          neg_count += AA_COUNT.Find( GetAATypes().ASP)->second;
        }
        if( AA_COUNT.Has( GetAATypes().GLU))
        {
          neg_count += AA_COUNT.Find( GetAATypes().GLU)->second;
        }

        storage::Row< double> &new_row( RESULT.InsertRow( "number_negative_aas"));
        new_row( 0) = neg_count;
        new_row( 1) = size_t( ( neg_count * 1000) / number_aas) / 10.0;
      }

      // pos charged residues
      {
        size_t pos_count( 0);
        if( AA_COUNT.Has( GetAATypes().ARG))
        {
          pos_count += AA_COUNT.Find( GetAATypes().ARG)->second;
        }
        if( AA_COUNT.Has( GetAATypes().LYS))
        {
          pos_count += AA_COUNT.Find( GetAATypes().LYS)->second;
        }

        storage::Row< double> &new_row( RESULT.InsertRow( "number_positive_aas"));
        new_row( 0) = pos_count;
        new_row( 1) = size_t( ( pos_count * 1000) / number_aas) / 10.0;
      }
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_ramachandran.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "biol/bcl_biol_environment_types.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score.h"
#include "util/bcl_util_enumerated.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! @brief get the default static instance of this class
    //! @return static default instance of this class
    const Ramachandran &Ramachandran::GetDefaultInstance()
    {
      // initialize a static const instance
      static const Ramachandran s_default_instance( GetDefaultHistogramFilename());

      // end
      return s_default_instance;
    }

    //! @brief get the default SSTypeHistogramFilename
    //! @return the default SSTypeHistogramFilename
    const std::string &Ramachandran::GetDefaultHistogramFilename()
    {
      // initialize static const string to hold the default histogram filename
      static const std::string s_default_sstype_histogram_filename( "phi_psi_angles_by_sstype.histogram2D");

      // end
      return s_default_sstype_histogram_filename;
    }

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Ramachandran::s_Instance
    (
      util::Enumerated< Ramachandran>::AddInstance( new Ramachandran())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Ramachandran::Ramachandran()
    {
    }

    //! @brief constructor from a AAType and SSType histogram filename
    //! @param SS_TYPE_HISTOGRAM_FILENAME filename for the phi/psi histogram according to SSTypes and AATypes
    Ramachandran::Ramachandran( const std::string &SS_TYPE_HISTOGRAM_FILENAME) :
      m_HistogramFilename( SS_TYPE_HISTOGRAM_FILENAME),
      m_AATypeMap(),
      m_SSTypeMap()
    {
      // initialize the members
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief Clone function
    //! @return pointer to new Ramachandran
    Ramachandran *Ramachandran::Clone() const
    {
      return new Ramachandran( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Ramachandran::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the HistogramFilename
    //! @return the HistogramFilename
    const std::string &Ramachandran::GetHistogramFilename() const
    {
      return m_HistogramFilename;
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &Ramachandran::GetAlias() const
    {
      static const std::string s_alias( "Ramachandran");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Ramachandran::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Distribution of phi/psi angles.");
      serializer.AddInitializer
      (
        "histogram filename",
        "path to the histogram file representing the distribution",
        io::Serialization::GetAgent( &m_HistogramFilename),
        GetDefaultHistogramFilename()
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return a random phi psi for the given AAType
    //! @param AA_TYPE AAType of interest
    //! @return pair of random phi and psi values for the given AAType
    storage::VectorND< 2, double> Ramachandran::GetRandomPhiPsi( const AAType &AA_TYPE) const
    {
      // find the corresponding member in the AAType map
      storage::Map< AAType, HistogramDistributionPair>::const_iterator map_itr( m_AATypeMap.Find( AA_TYPE));

      // assert that it is found
      BCL_Assert( map_itr != m_AATypeMap.End(), "no distribution stored for AAType " + AA_TYPE.GetName());

      // get a random value
      return
        map_itr->second.Second().DetermineRandomCase2D
        (
          map_itr->second.First().GetBoundariesX().First(),
          map_itr->second.First().GetBoundariesY().First(),
          map_itr->second.First().GetBinSizeXY().First(),
          map_itr->second.First().GetBinSizeXY().Second()
        );
    }

    //! @brief return a random phi psi for the given AAType and SSType
    //! @param AA_TYPE AAType of interest
    //! @param SS_TYPE SSType of interest
    //! @return pair of random phi and psi values for the given AAType and SSType
    storage::VectorND< 2, double> Ramachandran::GetRandomPhiPsi
    (
      const AAType &AA_TYPE,
      const SSType SS_TYPE
    ) const
    {
      return this->GetRandomPhiPsi( AA_TYPE, SS_TYPE, GetEnvironmentTypes().e_Solution);
    }

    //! @brief return a random phi psi for the given AAType and SSType
    //! @param AA_TYPE AAType of interest
    //! @param SS_TYPE SSType of interest
    //! @return pair of random phi and psi values for the given AAType and SSType
    storage::VectorND< 2, double> Ramachandran::GetRandomPhiPsi
    (
      const AAType &AA_TYPE,
      const SSType SS_TYPE,
      const EnvironmentType &ENV_TYPE
    ) const
    {
      const storage::Map< SSType, storage::Map< AAType, HistogramDistributionPair> > &map
      (
        ENV_TYPE == GetEnvironmentTypes().e_MembraneCore ? m_SSTypeMapMembrane : m_SSTypeMap
      );

      // find the corresponding member in the AAType map
      storage::Map< SSType, storage::Map< AAType, HistogramDistributionPair> >::const_iterator
        map_itr_a( map.Find( SS_TYPE));

      // assert that it is found
      BCL_Assert( map_itr_a != map.End(), "no distribution stored for SSType " + SS_TYPE.GetName());

      // now find the distribution for the corresponding aatype
      storage::Map< AAType, HistogramDistributionPair>::const_iterator map_itr_b( map_itr_a->second.Find( AA_TYPE));

      // assert that a probability distribution is stored for this amino acid type
      BCL_Assert( map_itr_b != map_itr_a->second.End(), "no distribution stored for amino acid " + AA_TYPE.GetName());

      // get a random value
      return
        map_itr_b->second.Second().DetermineRandomCase2D
        (
          map_itr_b->second.First().GetBoundariesX().First(),
          map_itr_b->second.First().GetBoundariesY().First(),
          map_itr_b->second.First().GetBinSizeXY().First(),
          map_itr_b->second.First().GetBinSizeXY().Second()
        );
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool Ramachandran::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      if( !command::CommandState::IsInStaticInitialization())
      {
        Initialize();
      }
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief function to read the ramachandran distributions and initialize the class members
    void Ramachandran::Initialize()
    {
      // reset
      m_SSTypeMap.Reset();
      m_AATypeMap.Reset();

      // initialize read
      io::IFStream read;

      // open read to the SSType PhiPsi histograms
      io::File::MustOpenIFStream( read, score::Score::AddHistogramPath( m_HistogramFilename));

      // iterate over the SSTypes
      for( size_t env_index( 0), n_env( 2); env_index < n_env; ++env_index)
      {
        std::string tmp;
        read >> tmp;
        for
        (
          SSTypes::const_iterator sstype_itr( GetSSTypes().Begin()), sstype_itr_end( GetSSTypes().COIL.GetIterator() + 1);
            sstype_itr != sstype_itr_end; ++sstype_itr
        )
        {
          // read in the enum
          SSType current_sstype;
          read >> current_sstype;

          // make sure it's correct
          BCL_Assert
          (
            *sstype_itr == current_sstype, "Wrong SSType " + sstype_itr->GetName() + " vs " + current_sstype.GetName()
          );

          // iterate over the valid AATypes
          for
          (
            AATypes::const_iterator
              aatype_itr( GetAATypes().Begin()),
              aatype_itr_end( GetAATypes().GetEnumIteratorFromIndex( AATypes::s_NumberStandardAATypes));
            aatype_itr != aatype_itr_end; ++aatype_itr
          )
          {
            // initialize pair to hold the histogram distribution pair
            HistogramDistributionPair hist_dist_pair;

            // read one letter code that represents the aatype
            std::string tmp;
            read >> tmp;

            // get the current amino acid type from the one letter code in "tmp"
            const AAType current_aatype( GetAATypes().AATypeFromOneLetterCode( tmp[ 0]));

            // assert that the aatypes in the histogram file are in the same order
            BCL_Assert
            (
              *aatype_itr == current_aatype,
              "unexpected aatype read from file! " + aatype_itr->GetName() + " != " +   current_aatype->GetName()
            );

            // read Histogram2D from stream
            read >> hist_dist_pair.First();

            if( !env_index)
            {
              // add to sstype independent map
              m_AATypeMap[ current_aatype].First().Combine( hist_dist_pair.First());
            }

            // create a Histogram2DDistribution from the histogram and store it
            hist_dist_pair.Second() = random::Histogram2DDistribution( hist_dist_pair.First());

            // insert into the member variable m_SSTypeMap
            ( env_index ? m_SSTypeMapMembrane : m_SSTypeMap)[ current_sstype][ current_aatype] = hist_dist_pair;
          }
        }
      }

      // clear the stream
      io::File::CloseClearFStream( read);

      // iterate over sstype independent map
      for
      (
        storage::Map< AAType, HistogramDistributionPair>::iterator
          itr( m_AATypeMap.Begin()), itr_end( m_AATypeMap.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->second.Second() = random::Histogram2DDistribution( itr->second.First());
      }
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_rotamer.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Rotamer::s_Instance
    (
      GetObjectInstances().AddInstance( new Rotamer())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Rotamer::Rotamer() :
      m_ChiAngles()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Rotamer
    Rotamer *Rotamer::Clone() const
    {
      return new Rotamer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Rotamer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief determines if rotamer contains any chi angles or not
    //! @return boolean - true if the rotamer contains no chi angles - false otherwise
    bool Rotamer::IsEmpty() const
    {
      return m_ChiAngles.IsEmpty();
    }

    //! @brief gives iterator to first chi in rotamer
    //! @return iterator to first chi in rotamer
    Rotamer::const_iterator Rotamer::Begin() const
    {
      return m_ChiAngles.Begin();
    }

    //! @brief gives iterator to end
    //! @return iterator to end
    Rotamer::const_iterator Rotamer::End() const
    {
      return m_ChiAngles.End();
    }

    //! @brief gives the number of chi angles in the rotamer
    //! @return the number of chi angles in the rotamer
    size_t Rotamer::GetSize() const
    {
      return m_ChiAngles.GetSize();
    }

    //! @brief gives the set of chis contained in this rotamer
    //! @return the set of chis contained in this rotamer
    storage::Set< ChiAngle::ChiEnum> Rotamer::GetChis() const
    {
      // to hold the chis contained in this rotamer
      storage::Set< ChiAngle::ChiEnum> chis;

      // iterate through the chi angles in order to fill the set of chis
      for
      (
        Rotamer::const_iterator chi_itr( Begin()), chi_itr_end( End());
        chi_itr != chi_itr_end;
        ++chi_itr
      )
      {
        chis.Insert( chi_itr->GetChi());
      }

      return chis;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief inserts a chi angle into the rotamer
    //! @param CHI_ANGLE the chi angle to add to the rotamer
    //! @return pair iterator to position of inserted element and bool indicating success or not
    std::pair< Rotamer::const_iterator, bool> Rotamer::Insert( const ChiAngle &CHI_ANGLE)
    {
      return m_ChiAngles.Insert( CHI_ANGLE);
    }

    //! @brief gives the angle value of the desired chi
    //! @param CHI the chi whose angle value is desired
    //! @param ANGLE_UNIT the unit the angle should be given in
    //! @return double which is the angle of the desired chi in desird units - undefined if chi does not exist
    double Rotamer::GetAngle( const ChiAngle::Chi &CHI, const math::Angle::Unit &ANGLE_UNIT) const
    {
      // try to find the chi angle with CHI
      const_iterator chi_itr( m_ChiAngles.Find( ChiAngle( CHI)));

      // true if chi could not be found
      if( chi_itr == m_ChiAngles.End())
      {
        return util::GetUndefinedDouble();
      }

      return chi_itr->GetAngle( ANGLE_UNIT);
    }

    //! @brief determines which chi angles have the same value between this and a given rotamer
    //!        in order for a chi to match, all previous chi also must match
    //! @param ROTAMER the other Rotamer whose chi angles will be compared to this
    //! @param ANGLE_UNIT the unit the angle tolerance is provided in
    //! @return set with all Chi that match dependent on the previous chi also being equivalent
    storage::Set< ChiAngle::ChiEnum> Rotamer::ChiMatchDependent
    (
      const Rotamer &ROTAMER, const math::Angle::Unit &ANGLE_UNIT, const double TOLERANCE
    ) const
    {
      storage::Set< ChiAngle::ChiEnum> matching_chi;

      for
      (
        Rotamer::const_iterator this_itr( Begin()), this_itr_end( End()),
          other_itr( ROTAMER.Begin()), other_itr_end( ROTAMER.End());
        this_itr != this_itr_end && other_itr != other_itr_end;
        ++this_itr, ++other_itr
      )
      {
        const double chi_angle_diff( this_itr->CalculateAngleDifference( *other_itr, ANGLE_UNIT));

        BCL_MessageDbg( "chi_angle_diff is " + util::Format()( chi_angle_diff));
        if( chi_angle_diff < TOLERANCE && util::IsDefined( chi_angle_diff))
        {
          matching_chi.Insert( this_itr->GetChi());
        }
        else
        {
          break;
        }
      }

      return matching_chi;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Rotamer::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ChiAngles, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Rotamer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ChiAngles, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief reads stream formatted in way easier for user to create
    //!        format is space separated on each line
    //!        <class identifier>
    //!        <chi enum> <angle value> <angle unit>
    //!        <chi enum> <angle value> <angle unit>
    //!        continued for as many chi angles as needed for this rotamer.
    //!        An example is
    //!        bcl::biol::Rotamer
    //!        e_Two 90 degree
    //!        e_Three 1.1 radian
    //! @return istream the rotamer was read from
    std::istream &Rotamer::ReadSimple( std::istream &ISTREAM)
    {
      BCL_MessageDbg( "reading in identifier");
      util::ObjectInterface::ReadIdentifier( ISTREAM);
      BCL_MessageDbg( "done reading in identifier");
      std::string line;
      while
      (
        !ISTREAM.eof() && line != GetClassIdentifier()
      )
      {
        std::getline( ISTREAM, line);
        util::TrimString( line);
        std::stringstream read( line);
        ChiAngle current_angle;
        if( !line.empty() && line != GetClassIdentifier())
        {
          current_angle.ReadSimple( read);
          BCL_MessageDbg( "current line is |" + line + "| and GetClassIdentifier is " + GetClassIdentifier());
          BCL_Assert
          (
            m_ChiAngles.Insert( current_angle).second, "could not insert chi angle " + util::Format()( current_angle) +
            "\ninto rotamer\n" + util::Format()( *this)
          );
        }
      }
        BCL_MessageDbg( "out of while loop current line is |" + line);

      return ISTREAM;
    }

    //! @brief gives description of this in format as read by ReadSimple
    //! @return std::string gives description of this in ReadSimple format
    std::string Rotamer::WriteSimple( const math::Angle::Unit &ANGLE_UNIT) const
    {
      std::string identification( "\n" + GetClassIdentifier());
      for
      (
        Rotamer::const_iterator itr( Begin()), itr_end( End());
        itr != itr_end;
        ++itr
      )
      {
        identification +=
        (
          "\n" + itr->WriteSimple( ANGLE_UNIT)
        );
      }

      return identification;
    }

    //! @brief binary functor helper struct for determing if one chi angle is less than another
    //!        sorts according to the ChiEnum within the chi angles
    //! @param CHI_ANGLE_A first chi angle object
    //! @param CHI_ANGLE_B second chi angle object
    //! @return bool true if the chi of CHI_ANGLE_A is less than CHI_ANGLE_B - false otherwise
    bool Rotamer::ChiAngleLessThan::operator()( const ChiAngle &CHI_ANGLE_A, const ChiAngle &CHI_ANGLE_B) const
    {
      return CHI_ANGLE_A.GetChi() < CHI_ANGLE_B.GetChi();
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_rotamer_library.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> RotamerLibrary::s_Instance
    (
      GetObjectInstances().AddInstance( new RotamerLibrary())
    );

    //! map of loop libraries
    storage::HashMap< std::string, util::ShPtr< RotamerLibrary> > RotamerLibrary::s_RotamerLibraries =
      storage::HashMap< std::string, util::ShPtr< RotamerLibrary> >();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RotamerLibrary::RotamerLibrary() :
      m_LibraryFileName()
    {
    }

    //! @brief construct from members
    //! @param LIBRARY_FILE_NAME path to the rotamer library
    RotamerLibrary::RotamerLibrary( const std::string &LIBRARY_FILE_NAME) :
      m_LibraryFileName( LIBRARY_FILE_NAME)
    {
      // read in the rotamer library at the given path
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief returns a loop library from the given file name
    //! @param LIBRARY_FILE_NAME path to the rotamer library
    util::ShPtr< RotamerLibrary> RotamerLibrary::CreateRotamerLibrary( const std::string &LIBRARY_FILE_NAME)
    {
      if( command::CommandState::GetGlobalCommandState().IsInStaticInitialization())
      {
        return util::ShPtr< RotamerLibrary>();
      }

      // compute hash key for the given parameters
      const std::string key( LIBRARY_FILE_NAME);

      // if a library with this key does not exist, create one
      if( s_RotamerLibraries.Find( key) == s_RotamerLibraries.End())
      {
        util::ShPtr< RotamerLibrary> sp_library( new RotamerLibrary( LIBRARY_FILE_NAME));
        s_RotamerLibraries[ key] = sp_library;
        return sp_library;
      }
      return s_RotamerLibraries[ key];
    }

    //! @brief copy constructor
    //! @return pointer to a new RotamerLibrary
    RotamerLibrary *RotamerLibrary::Clone() const
    {
      return new RotamerLibrary( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &RotamerLibrary::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &RotamerLibrary::GetAlias() const
    {
      static const std::string s_alias( "RotamerLibrary");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RotamerLibrary::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Library containing rotamers and corresponding probabilities for amino acid side chains.");
      serializer.AddInitializer
      (
        "file path",
        "path to the rotamer library file",
        io::Serialization::GetAgent( &m_LibraryFileName)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool RotamerLibrary::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      if( !command::CommandState::IsInStaticInitialization())
      {
        // read the rotamer library file
        io::IFStream lib_file;
        io::File::MustOpenIFStream( lib_file, m_LibraryFileName);
        ReadLibrary( lib_file);
        io::File::CloseClearFStream( lib_file);
      }

      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief read a loop template library from the given input stream
    //! @param ISTREAM input stream from which to read the loop template library
    //! @param ANGLE_BIN_WIDTH bin width for chi angles
    void RotamerLibrary::ReadLibrary( std::istream &ISTREAM, double ANGLE_BIN_WIDTH)
    {
      // read in all lines in the library
      storage::Vector< std::string> lines( util::StringLineListFromIStream( ISTREAM));

      // iterate over all lines and add each valid rotamer to the library
      storage::HashMap< std::string, storage::HashMap< size_t, size_t> > rotamer_probability;
      for( auto line_itr( lines.Begin()), line_itr_end( lines.End()); line_itr != line_itr_end; ++line_itr)
      {
        // skip empty lines or comment lines indicated by a leading # or !
        if( line_itr->empty() || ( *line_itr)[ 0] == '!' || ( *line_itr)[ 0] == '#')
        {
          continue;
        }

        // read in the relevant data for this rotamer
        AAType aa_type;                     // type of this amino acid
        size_t count;                             // number of side chains
        storage::Vector< size_t> rotamers( 4, 0); // rotamer configuration
        storage::Vector< std::string> split_line( util::SplitString( *line_itr, " \t"));
        std::string aa_type_tmp( *split_line[ 0]);
        util::TryConvertFromString( count, *split_line[ 3], util::GetLogger());
        for( size_t rotamer_angles( 0); rotamer_angles < 4; ++rotamer_angles)
        {
          util::TryConvertFromString( rotamers( rotamer_angles), *split_line[ rotamer_angles + 4], util::GetLogger());
        }

        // in this implementation the backbone dihedral angles are ignored and the rotamer conformations are aggregated over
        // all backbone conformations
        size_t key( 0), radix( 1);
        for( size_t i( 0); i < 4; ++i, radix *= 10)
        {
          key += radix * rotamers( i);
        }
        rotamer_probability[ aa_type_tmp][ key] += count;
      }

    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_sasa_data.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "math/bcl_math_cubic_spline.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasaData::s_Instance
    (
      GetObjectInstances().AddInstance( new SasaData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasaData::SasaData() :
      m_SasaData()
    {
    }

    //! @brief constructor from given input data
    SasaData::SasaData( const storage::Vector< SasaPoint> &INIT_DATA) :
      m_SasaData( INIT_DATA)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SasaData
    SasaData *SasaData::Clone() const
    {
      return new SasaData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasaData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief pushback function to add object to Dataset vector
    //! @param DATAPOINT_OBJECT DataPoint values Q, I, and Error
    void SasaData::PushBack( const SasaPoint &DATAPOINT_OBJECT)
    {
      m_SasaData.PushBack( DATAPOINT_OBJECT);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief preallocate memory
    //! @param SIZE size to preallocate
    void SasaData::AllocateMemory( const size_t &SIZE)
    {
      m_SasaData.AllocateMemory( SIZE);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads in the member data from a formatted file containing 3 columns:
    //! @param ISTREAM input stream
    //! @param FORMAT the file format to use for reading
    //! @return istream which was read from
    std::istream &SasaData::ReadFromDataFile( std::istream &ISTREAM)
    {
      // build a string to hold the line information
      std::string read_line;

      // skip the first line
      std::getline( ISTREAM, read_line);

      // while the end of the file is not reached
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        // push back the data
        m_SasaData.PushBack
        (
          SasaPoint
            (
              util::ConvertStringToNumericalValue< double>( split_line( 0)),    // Atom Number
              util::ConvertStringToNumericalValue< double>( split_line( 1)),    // Solvent Excluded Surface
              util::ConvertStringToNumericalValue< double>( split_line( 2))     // Solvent Accessible Surface
            )
         );
        }

      size_t size( m_SasaData.GetSize());

      BCL_Assert( size != 0, "The number of Sasa values are: " + util::Format()( size));

      // return the stream
      return ISTREAM;
    }

    //! @brief writes out the member data from a formatted file containing 3 columns
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &SasaData::WriteToDataFile( std::ostream &OSTREAM) const
    {
      // initialize format
      const util::Format format;

      // iterate over m_SasaData
      for
      (
        storage::Vector< SasaPoint>::const_iterator itr( m_SasaData.Begin()), itr_end( m_SasaData.End());
         itr != itr_end; ++itr
      )
      {
        // write the data
        OSTREAM << format( itr->GetAtomNumber()) << '\t'
                << format( itr->GetSolventExcludedSurface()) << '\t'
                << format( itr->GetSolventAccessibleSurface()) << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasaData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SasaData, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasaData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SasaData, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_sasa_point.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasaPoint::s_Instance
    (
      GetObjectInstances().AddInstance( new SasaPoint())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasaPoint::SasaPoint() :
      m_AtomNumber( util::GetUndefinedSize_t()),
      m_SolventExcludedSurface( util::GetUndefinedDouble()),
      m_SolventAccessibleSurface( util::GetUndefinedDouble())
    {
    }

    //! @brief constructor from given input data
    //! @param ATOMNUM - The atom number in the PDB
    //! @param SES     - The Solvent Excluded Surface
    //! @param SAS     - The Solvent Accessible Surface
    SasaPoint::SasaPoint
    (
      const size_t ATOMNUM,
      const double SES,
      const double SAS
    ) :
       m_AtomNumber( ATOMNUM),
       m_SolventExcludedSurface( SES),
       m_SolventAccessibleSurface( SAS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SasaPoint
    SasaPoint *SasaPoint::Clone() const
    {
      return new SasaPoint( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasaPoint::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasaPoint::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AtomNumber, ISTREAM);
      io::Serialize::Read( m_SolventExcludedSurface, ISTREAM);
      io::Serialize::Read( m_SolventAccessibleSurface, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasaPoint::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AtomNumber, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SolventExcludedSurface, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SolventAccessibleSurface, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_ss_type_data.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "math/bcl_math_rotation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSTypeData::s_Instance( GetObjectInstances().AddInstance( new SSTypeData()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct undefined secondary structure type
    SSTypeData::SSTypeData() :
      m_OneLetterCode( ' '),
      m_IsStructured( false),
      m_RadialExtent( util::GetUndefined< double>()),
      m_AnglePerTurn( util::GetUndefined< double>()),
      m_RiseInZPerResidue( util::GetUndefined< double>()),
      m_IdealPhi( util::GetUndefined< double>()),
      m_IdealPsi( util::GetUndefined< double>()),
      m_FragmentLength( util::GetUndefined< size_t>()),
      m_ContactWindowRadius( util::GetUndefined< size_t>()),
      m_TransformationMatrixForResidues(),
      m_ThreeStatePrediction( util::GetUndefined< double>()),
      m_PhiRange(),
      m_PsiRange()
    {
    }

    //! @brief construct secondary structure type from provided data
    //! @param ONE_LETTER_CODE one letter description
    //! @param IS_STRUCTURED Whether it has a defined geometry (only helix and strand has it)
    //! @param RADIAL_EXTENT radial extensions in angstroms
    //! @param ANGLE_PER_TURN angular turn per residue
    //! @param RISE_IN_Z_PER_RESIDUE rise in z axis per residue
    //! @param IDEAL_PHI ideal phi angle
    //! @param IDEAL_PSI ideal psi angle
    //! @param FRAGMENT_LENGTH Length of fragments created from this SSType that describe the curvature
    //! @param CONTACT_WINDOW_RADIUS number of neighbor residues in contact prediction descriptors
    //! @param THREE_STATE_PREDICTION Vector that contains three state prediction
    //! @param BACKBONE_PHI_RANGE permissible range of backbone phi angle for this ss type
    //! @param BACKBONE_PSI_RANGE permissible range of backbone psi angle for this ss type
    SSTypeData::SSTypeData
    (
      const char ONE_LETTER_CODE,
      const bool IS_STRUCTURED,
      const double RADIAL_EXTENT,
      const double ANGLE_PER_TURN,
      const double RISE_IN_Z_PER_RESIDUE,
      const double IDEAL_PHI,
      const double IDEAL_PSI,
      const size_t FRAGMENT_LENGTH,
      const size_t CONTACT_WINDOW_RADIUS,
      const linal::Vector3D &THREE_STATE_PREDICTION,
      const math::Range< double> &BACKBONE_PHI_RANGE,
      const math::Range< double> &BACKBONE_PSI_RANGE
    ) :
      m_OneLetterCode( ONE_LETTER_CODE),
      m_IsStructured( IS_STRUCTURED),
      m_RadialExtent( RADIAL_EXTENT),
      m_AnglePerTurn( ANGLE_PER_TURN),
      m_RiseInZPerResidue( RISE_IN_Z_PER_RESIDUE),
      m_IdealPhi( IDEAL_PHI),
      m_IdealPsi( IDEAL_PSI),
      m_FragmentLength( FRAGMENT_LENGTH),
      m_ContactWindowRadius( CONTACT_WINDOW_RADIUS),
      m_TransformationMatrixForResidues
      (
        math::TransformationMatrix3D
        (
          linal::Vector3D( 0.0, 0.0, RISE_IN_Z_PER_RESIDUE)
        )( math::RotationMatrix3D( coord::GetAxes().e_Z, ANGLE_PER_TURN))
      ),
      m_ThreeStatePrediction( THREE_STATE_PREDICTION),
      m_PhiRange( BACKBONE_PHI_RANGE),
      m_PsiRange( BACKBONE_PSI_RANGE)
    {
    }

    //! @brief virtual copy constructor
    SSTypeData *SSTypeData::Clone() const
    {
      return new SSTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return one letter code
    //! @return one letter code
    const char SSTypeData::GetOneLetterCode() const
    {
      return m_OneLetterCode;
    }

    //! @brief return if it is's structured
    //! @return whether it's structured
    bool SSTypeData::IsStructured() const
    {
      return m_IsStructured;
    }

    //! @brief return radial extensions in angstroms
    //! @return radial extensions in angstroms
    const double SSTypeData::GetRadialExtent() const
    {
      return m_RadialExtent;
    }

    //! @brief return angle per turn
    //! @return angle per turn
    const double SSTypeData::GetAnglePerTurn() const
    {
      return m_AnglePerTurn;
    }

    //! @brief return rise in z axis per residue
    //! @return rise in z axis per residue
    const double SSTypeData::GetRiseInZPerResidue() const
    {
      return m_RiseInZPerResidue;
    }

    //! @brief return ideal phi angle
    //! @return ideal phi angle
    const double SSTypeData::GetIdealPhi() const
    {
      return m_IdealPhi;
    }

    //! @brief return ideal psi angle
    //! @return ideal psi angle
    const double SSTypeData::GetIdealPsi() const
    {
      return m_IdealPsi;
    }

    //! @brief return Length of fragments created from this SSType that describe the curvature
    //! @return Length of fragments created from this SSType that describe the curvature
    const size_t SSTypeData::GetFragmentLength() const
    {
      return m_FragmentLength;
    }

    //! @brief return number of neighboring residues on each side to be used in contact prediction descriptors
    //! @return number of neighboring residues on each side to be used in contact prediction descriptors
    const size_t SSTypeData::GetContactWindowRadius() const
    {
      return m_ContactWindowRadius;
    }

    //! @brief return total number of neighbor residues to be used in contact prediction descriptors
    //! @return total number of neighbor residues to be used in contact prediction descriptors
    const size_t SSTypeData::GetContactWindowLength() const
    {
      return 2 * m_ContactWindowRadius + 1;
    }

    //! @brief return transformation matrix for residues for this SSType
    //! @return transformation matrix for residues for this SSType
    const math::TransformationMatrix3D &SSTypeData::GetTransformationMatrixForResidues() const
    {
      return m_TransformationMatrixForResidues;
    }

    //! @brief returns three state prediction vector with this SSType with a prediction of 1.0
    //! @return three state prediction vector with this SSType with a prediction of 1.0
    const linal::Vector3D &SSTypeData::GetThreeStatePrediction() const
    {
      return m_ThreeStatePrediction;
    }

    //! @brief returns backbone phi range for this SSType
    //! @return range of backbone phi angles that are permissible
    const math::Range< double> SSTypeData::GetBackbonePhiRange() const
    {
      return m_PhiRange;
    }

    //! @brief returns backbone psi range for this SSType
    //! @return range of backbone psi angles that are permissible
    const math::Range< double> SSTypeData::GetBackbonePsiRange() const
    {
      return m_PsiRange;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_OneLetterCode,        ISTREAM);
      io::Serialize::Read( m_IsStructured,         ISTREAM);
      io::Serialize::Read( m_RadialExtent,         ISTREAM);
      io::Serialize::Read( m_AnglePerTurn,         ISTREAM);
      io::Serialize::Read( m_RiseInZPerResidue,    ISTREAM);
      io::Serialize::Read( m_IdealPhi,             ISTREAM);
      io::Serialize::Read( m_IdealPsi,             ISTREAM);
      io::Serialize::Read( m_FragmentLength,       ISTREAM);
      io::Serialize::Read( m_ContactWindowRadius,  ISTREAM);
      io::Serialize::Read( m_ThreeStatePrediction, ISTREAM);
      io::Serialize::Read( m_PhiRange,             ISTREAM);
      io::Serialize::Read( m_PsiRange,             ISTREAM);

      // initialize transformation
      m_TransformationMatrixForResidues = math::TransformationMatrix3D( linal::Vector3D( 0.0, 0.0, m_RiseInZPerResidue));
      m_TransformationMatrixForResidues( math::RotationMatrix3D( coord::GetAxes().e_Z, m_AnglePerTurn));

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &SSTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_OneLetterCode,        OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IsStructured,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RadialExtent,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AnglePerTurn,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RiseInZPerResidue,    OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IdealPhi,             OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IdealPsi,             OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FragmentLength,       OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ContactWindowRadius,  OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ThreeStatePrediction, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PhiRange,             OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PsiRange,             OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
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
#include "biol/bcl_biol_ss_types.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_angle.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief constructor all SSTypes
    SSTypes::SSTypes() :
//                                                               IsStructured                  Radial Extent                   AnglePerTurn                  RiseInZPerRes                            Phi                            Psi                     FragLength            ContactWindowRadius            ThreeStatePrediction      backbone phi range      backbone psi range
      HELIX(              AddEnum( "HELIX"           , SSTypeData( 'H',  true,                         4.240,  math::Angle::Radian( -100.0),                       1.50247,                     -0.996767,                      -0.80718,                             5,                             4, linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>( math::Angle::Radian( -135.0), math::Angle::Radian( -25.0)), math::Range< double>( math::Angle::Radian( -70.0), math::Angle::Radian(  20.0))))),
      STRAND(             AddEnum( "STRAND"          , SSTypeData( 'E',  true,                         3.275,  math::Angle::Radian( -180.0),                         3.425,                      -2.39835,                       2.35502,                             3,                             2, linal::Vector3D( 0.0, 1.0, 0.0), math::Range< double>( math::Angle::Radian( -180.0), math::Angle::Radian( -35.0)), math::Range< double>( math::Angle::Radian(  25.0), math::Angle::Radian( 180.0))))),
      COIL(               AddEnum( "COIL"            , SSTypeData( 'C', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 0.0, 0.0, 1.0), math::Range< double>(), math::Range< double>()))),
      e_HelixRightOmega(  AddEnum( "HelixRightOmega" , SSTypeData( '2', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix right handed omega
      e_HelixRightPi(     AddEnum( "HelixRightPi"    , SSTypeData( '3', false, util::GetUndefined< double>(),  math::Angle::Radian(  -87.0),                          1.15, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix right handed pi
      e_HelixRightGamma(  AddEnum( "HelixRightGamma" , SSTypeData( '4', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix right handed gamma
      e_HelixRight310(    AddEnum( "HelixRight310"   , SSTypeData( '5', false, util::GetUndefined< double>(),  math::Angle::Radian( -120.0),                           2.0, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix right handed 3-10
      e_HelixLeftAlpha(   AddEnum( "HelixLeftAlpha"  , SSTypeData( '6', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>( math::Angle::Radian(   45.0), math::Angle::Radian(  65.0)), math::Range< double>( math::Angle::Radian(  15.0), math::Angle::Radian( 100.0))))), // helix left handed alpha
      e_HelixLeftOmega(   AddEnum( "HelixLeftOmega"  , SSTypeData( '7', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix left handed omega
      e_HelixLeftGamma(   AddEnum( "HelixLeftGamma"  , SSTypeData( '8', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix left handed gamma
      e_Helix27Ribbon(    AddEnum( "Helix27Ribbon"   , SSTypeData( '9', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix 2 - 7 ribbon
      e_HelixPolyProline( AddEnum( "HelixPolyProline", SSTypeData( '0', false, util::GetUndefined< double>(),  math::Angle::Radian(  120.0),                           3.1, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>( math::Angle::Radian(  -80.0), math::Angle::Radian( -70.0)), math::Range< double>( math::Angle::Radian( 145.0), math::Angle::Radian( 155.0)))))  // helix poly proline
    {
      m_HelixTypes = storage::Set< SSType>( Begin() + e_HelixRightOmega, Begin() + e_HelixPolyProline + 1);
      m_HelixTypes.Insert( HELIX);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the reduced (3-state) types
    const storage::Vector< SSType> &SSTypes::GetReducedTypes() const
    {
      static storage::Vector< SSType> s_reduced( storage::Vector< SSType>::Create( HELIX, STRAND, COIL));
      return s_reduced;
    }

    //! @brief Get secondary structure type from provided ONE_LETTER_CODE
    //! @param ONE_LETTER_CODE one letter code for SSType
    //! @return secondary structure type from provided ONE_LETTER_CODE
    const SSType &SSTypes::SSTypeFromOneLetterCode( const char ONE_LETTER_CODE) const
    {
      // iterate over all ss types
      for
      (
        const_iterator type_itr( Begin()), type_itr_end( End());
        type_itr != type_itr_end;
        ++type_itr
      )
      {
        // if given ONE_LETTER_CODE matches this one letter code
        if( ( *type_itr)->GetOneLetterCode() == ONE_LETTER_CODE)
        {
          // return
          return *type_itr;
        }
      }

      // if no match is found return undefined
      return GetSSTypes().e_Undefined;
    }

    //! @brief secondary structure from phi and psi angle
    //! @param PHI phi backbone angle
    //! @param PSI psi backbone angle
    //! @return most likely SSType - if non was found, COIL
    const SSType &SSTypes::SSTypeFromPhiPsi( const double PHI, const double PSI) const
    {
      for( const_iterator sstype_itr( Begin()), sstype_itr_end( End()); sstype_itr != sstype_itr_end; ++sstype_itr)
      {
        if( ( *sstype_itr)->GetBackbonePhiRange().IsWithin( PHI) && ( *sstype_itr)->GetBackbonePsiRange().IsWithin( PSI))
        {
          return *sstype_itr;
        }
      }

      return COIL;
    }

    //! @brief converts given pdb helix class into SSType
    //! @see http://www.wwpdb.org/documentation/format32/sect5.html#HELIX
    //! @param HELIX_CLASS as found in HelixClass entry in pdb HELIX line 1-10
    //! @return helix SSType for given class, e_Undefined if HELIX_CLASS is not in [1,10] range
    const SSType &SSTypes::SSTypeFromPDBHelixClass( const size_t HELIX_CLASS) const
    {
      // from: http://www.wwpdb.org/documentation/format32/sect5.html#HELIX
      // TYPE OF  HELIX                     (COLUMNS 39 - 40)
      // --------------------------------------------------------------
      // Right-handed alpha (default)                1
      // Right-handed omega                          2
      // Right-handed pi                             3
      // Right-handed gamma                          4
      // Right-handed 3 - 10                         5
      // Left-handed alpha                           6
      // Left-handed omega                           7
      // Left-handed gamma                           8
      // 2 - 7 ribbon/helix                          9
      // Polyproline                                10
      static const math::Range< size_t> s_valid_classes( 1, 10);

      if( !s_valid_classes.IsWithin( HELIX_CLASS))
      {
        return HELIX;
//        return e_Undefined;
      }

      // right handed alpha helix
      if( HELIX_CLASS == 1)
      {
        return HELIX;
      }

      // remaining helices
      return *( e_HelixRightOmega.GetIterator() + ( HELIX_CLASS - s_AlphaOmegaHelixOffset));
    }

    //! @brief converts a SSType to a pdb helix class
    //! @param SSTYPE the sstype to convert to helix class
    //! @return pdb helix class [1-10]; if SSTYPE is not a valid helix class, undefined
    size_t SSTypes::PDBHelixClassFromSSType( const SSType &SSTYPE) const
    {
      // alpha helix is 1
      if( SSTYPE == HELIX)
      {
        return 1;
      }

      // strand coil undefined are undefined
      if( !SSTYPE.IsDefined() || SSTYPE <= COIL)
      {
        return util::GetUndefined< size_t>();
      }

      // rest is directly related to ss type
      const size_t helix_class( SSTYPE.GetIndex() + 1 - s_AlphaOmegaHelixOffset);
      if( helix_class > 10)
      {
        return util::GetUndefined< size_t>();
      }

      return helix_class;
    }

    //! @brief set of all helix types
    //! @return Set containing all helix types
    const storage::Set< SSType> &SSTypes::GetHelixTypes() const
    {
      return m_HelixTypes;
    }

    //! @brief test if two sses are of similar type (alpha helix to 3-10 helix are considered similar, but strand would not be)
    //! @param SS_TYPE_LHS left hand side ss type
    //! @param SS_TYPE_RHS right hand side ss type
    //! @return true if the two types are considered similar
    bool SSTypes::AreSimilar( const SSType &SS_TYPE_LHS, const SSType &SS_TYPE_RHS) const
    {
      // identical type is similar
      if( SS_TYPE_LHS == SS_TYPE_RHS)
      {
        return true;
      }

      // consider right helices as similar
      if( SS_TYPE_LHS == HELIX)
      {
        if( SS_TYPE_RHS >= e_HelixRightOmega && SS_TYPE_RHS <= e_HelixRight310)
        {
          return true;
        }
      }
      if( SS_TYPE_RHS == HELIX)
      {
        if( SS_TYPE_LHS >= e_HelixRightOmega && SS_TYPE_LHS <= e_HelixRight310)
        {
          return true;
        }
      }

      return false;
    }

    //! @brief on access function for all SSTypes
    //! @return SSTypes enums
    const SSTypes &GetSSTypes()
    {
      return SSTypes::GetEnums();
    }

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< biol::SSTypeData, biol::SSTypes>;

  } // namespace util
} // namespace bcl
