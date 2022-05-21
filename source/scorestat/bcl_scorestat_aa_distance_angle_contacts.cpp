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
// include header for this class
#include "scorestat/bcl_scorestat_aa_distance_angle_contacts.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_atom_sigma_charge.h"
#include "descriptor/bcl_descriptor_coulombic_force.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_2d.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_histogram_3d.h"
#include "math/bcl_math_running_average.h"
#include "score/bcl_score_aa_pair_sidechain_interaction.h"
#include "util/bcl_util_string_functions.h"
#include "util/bcl_util_wrapper.h"
namespace bcl
{
  namespace scorestat
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AADistanceAngleContacts::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AADistanceAngleContacts())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AADistanceAngleContacts::AADistanceAngleContacts() :
      m_BinSize( 1.0),
      m_ChainIds( ""),
      m_AADistSeqExcl( 2),
      m_VdwRadiusStats( false),
      m_MinCounts( 10)
    {
    }

    //! @brief virtual copy constructor
    AADistanceAngleContacts *AADistanceAngleContacts::Clone() const
    {
      return new AADistanceAngleContacts( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AADistanceAngleContacts::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AADistanceAngleContacts::GetOutFilePostfix() const
    {
      static const std::string s_name( "contact_probability");
      return s_name;
    }

    //! @brief returns the binsize for the histogram
    //! @return the binsize for the histogram
    const double &AADistanceAngleContacts::GetBinSize() const
    {
      return m_BinSize;
    }

    //! @brief returns chain id
    //! @return chain id
    const std::string &AADistanceAngleContacts::GetChainId() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AADistanceAngleContacts::GetAlias() const
    {
      static const std::string s_name( "AADistanceAngleContacts");
      return s_name;
    }

  //////////////
  // operator //
  //////////////

    namespace
    {
      //! @brief helper function to get the vdw radii of one of the strings used for vdw computation below
      double GetVdwRadiiOfType( const std::string &A)
      {
        return A[ 0] == 'O' ? 0.7 // backbone oxygen; real vdw is 1.2, but vdw excludes hydrogen bonds, which backbone O almost always makes
               : (
                   A[ 0] == 'N' ? 1.34 // backbone nitrogen
                   : (
                       A[ 0] == 'G' && A[ 2] == 'Y'
                       ? 0.50 // Glycine HA. Often makes hydrogen bonds, so this radius is deliberately close to the covalent radius of H
                       : 1.38 // C/CA/CB
                     )
                 );
      }
    }

    //! @brief add statistics about the given protein and chains to the internal table and histograms
    //! @param ENSEMBLE the protein ensemble the statistics of which to be computed
    std::string AADistanceAngleContacts::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // initializes maps for storing histograms of distances between all types of amino acid pairs
      std::vector
      <
        std::vector
        <
          std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram3D> >
        >
      > distance_map_background
        (
          3,
          std::vector< std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram3D> > >( 3)
        );
      std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram3D> > distance_map_contacts;
      std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram3D> > rand_bg;

      storage::Vector< storage::Vector< math::RunningAverage< double> > > fraction_beyond_12A
      (
        20,
        storage::Vector< math::RunningAverage< double> >( 20)
      );

      storage::Map< std::string, storage::Vector< double> > vdw_map;

      // initialize maps for storing aa distances
      // iterate over all aa types to set the first aa
      for
      (
        biol::AATypes::const_iterator aa_1_itr( biol::GetAATypes().Begin()), aa_itr_end( biol::GetAATypes().End());
          aa_1_itr != aa_itr_end;
        ++aa_1_itr
      )
      {
        // iterate over all aa types to set the second aa
        for( biol::AATypes::const_iterator aa_2_itr( aa_1_itr); aa_2_itr != aa_itr_end; ++aa_2_itr)
        {
          // skip non-natural amino acids
          if( !( *aa_1_itr)->IsNaturalAminoAcid() || !( *aa_2_itr)->IsNaturalAminoAcid())
          {
            continue;
          }

          std::pair< biol::AAType, biol::AAType> types( *aa_1_itr, *aa_2_itr);
          distance_map_contacts[ types] =
            util::ShPtr< math::Histogram3D>
            (
              new math::Histogram3D
              (
                storage::VectorND< 3, double>( 2.0, -1.0, -1.0),
                storage::VectorND< 3, double>( m_BinSize, 2.0 / 3.0, 2.0 / 3.0),
                storage::VectorND< 3, size_t>( 10.0 / m_BinSize, 3, 3),
                -m_MinCounts
              )
            );

          util::ShPtr< math::Histogram3D> aa_rand_bg
          (
            rand_bg[ types] =
              util::ShPtr< math::Histogram3D>
              (
                new math::Histogram3D
                (
                  storage::VectorND< 3, double>( 2.0, -1.0, -1.0),
                  storage::VectorND< 3, double>( m_BinSize, 2.0 / 3.0, 2.0 / 3.0),
                  storage::VectorND< 3, size_t>( 10.0 / m_BinSize, 3, 3),
                  2.0
                )
              )
          );
          // background from taking random ca-cb vectors at a random, evenly distributed distance
          // and angle with the constraints that the other vectors be on different SSEs (so not too close
          // in XY or within clashing distance (too close overall distance)
          {
            linal::Vector3D cb_a( 0.00, 0, 0);
            const double a_length
            (
              chemistry::BondLengths::GetBondLength
              (
                chemistry::GetAtomTypes().C_TeTeTeTe, // C-Alpha atom type
                1,                                    // single bond
                ( *aa_1_itr)->GetFirstSidechainAtomType()->GetElementType() == chemistry::GetElementTypes().e_Carbon
                ? chemistry::GetAtomTypes().C_TeTeTeTe
                : chemistry::GetAtomTypes().H_S
              )
            );
            const double b_length
            (
              chemistry::BondLengths::GetBondLength
              (
                chemistry::GetAtomTypes().C_TeTeTeTe, // C-Alpha atom type
                1,                                    // single bond
                ( *aa_2_itr)->GetFirstSidechainAtomType()->GetElementType() == chemistry::GetElementTypes().e_Carbon
                ? chemistry::GetAtomTypes().C_TeTeTeTe
                : chemistry::GetAtomTypes().H_S
              )
            );
            // van-der waals radius of C-Alpha
            const double ca_vdw
            (
              chemistry::GetAtomTypes().C_TeTeTeTe->GetAtomTypeProperty( chemistry::AtomTypeData::e_VdWaalsRadiusCSD)
            );
            const double cb_a_vdw( ( *aa_1_itr)->GetVdwRadiusToOtherAA( ( *aa_1_itr)->GetFirstSidechainAtomType()));
            const double cb_b_vdw( ( *aa_2_itr)->GetVdwRadiusToOtherAA( ( *aa_2_itr)->GetFirstSidechainAtomType()));
            const double cb_a_ca_vdw( cb_a_vdw + ca_vdw);
            const double cb_b_ca_vdw( cb_b_vdw + ca_vdw);
            const double cb_cb_vdw( cb_b_vdw + cb_a_vdw);
            const double ca_ca_vdw( ca_vdw + ca_vdw);
            linal::Vector3D ca_a( a_length, 0, 0);
            const double half_typical_sse_length( 7.0);
            const double typical_sse_min_d_cb( 3.6);
            const double typical_sse_min_d_ca( 3.3);

            for( int i( 0); i < 1000000; ++i)
            {
              ca_a.Rotate( math::RotationMatrix3D().SetRand());
              linal::Vector3D cb_b( cb_a);
              cb_b.SetRandomTranslation( 12.0);
              if( linal::Distance( cb_b, ca_a) < cb_b_ca_vdw || linal::Distance( cb_b, cb_a) < cb_cb_vdw)
              {
                --i;
                continue;
              }
              linal::Vector3D ca_b( b_length, 0, 0);
              ca_b.Rotate( math::RotationMatrix3D().SetRand());
              ca_b += cb_b;
              if( linal::Distance( ca_b, ca_a) < ca_ca_vdw || linal::Distance( ca_b, cb_a) < cb_a_ca_vdw)
              {
                --i;
                continue;
              }
              linal::Vector2D ca_b_2d( ca_b.X(), ca_b.Y());
              if
              (
                ca_b_2d.Norm() < typical_sse_min_d_ca
                && ca_b.Z() < half_typical_sse_length
                && ca_b.Z() > -half_typical_sse_length
              )
              {
                --i;
                continue;
              }
              linal::Vector2D cb_b_2d( cb_b.X(), cb_b.Y());
              if
              (
                cb_b_2d.Norm() < typical_sse_min_d_cb
                && cb_b.Z() < half_typical_sse_length
                && cb_b.Z() > -half_typical_sse_length
              )
              {
                --i;
                continue;
              }
              const double dist( linal::Distance( cb_b, cb_a));
              double angle( linal::ProjAngleCosinus( linal::UnitVector( ca_a, cb_a), linal::UnitVector( ca_a, cb_b)));
              double angleb( linal::ProjAngleCosinus( linal::UnitVector( ca_b, cb_b), linal::UnitVector( ca_b, cb_a)));
              aa_rand_bg->PushBack( dist, angleb, angle, 1.0);
            }
          }
          for( size_t i( 0), n_sse_types( 3); i < n_sse_types; ++i)
          {
            for( size_t j( 0); j < n_sse_types; ++j)
            {
              distance_map_background[ i][ j][ std::pair< biol::AAType, biol::AAType>( *aa_1_itr, *aa_2_itr)] =
                util::ShPtr< math::Histogram3D>
                (
                  new math::Histogram3D
                  (
                    storage::VectorND< 3, double>( 2.0, -1.0, -1.0),
                    storage::VectorND< 3, double>( m_BinSize, 2.0 / 3.0, 2.0 / 3.0),
                    storage::VectorND< 3, size_t>( 10.0 / m_BinSize, 3, 3),
                    2.0
                  )
                );
            }
          }
        }
      }

      size_t number_pairs( 0);
      linal::Vector< size_t> number_aas( 20, size_t( 0));
      linal::Vector< size_t> number_interactions( 20, size_t( 0));

      // iterate through all models in the ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
          protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        util::ShPtr< assemble::ProteinModel> sp_current_model( *protein_model_itr);
        const assemble::ProteinModel &protein_model( *sp_current_model);

        // get current pdb name before all the loops start
        const util::ShPtr< util::Wrapper< std::string> > &model_name_ptr(
          protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");

        // instantiate sequences CaCb from pdb for each chain one
        BCL_MessageStd( "building sequences from pdb chains on " + model_name);
        util::SiPtrVector< const biol::AASequence> aa_sequences( protein_model.GetSequences());

        // number of chains
        BCL_MessageStd( "pdb has " + util::Format()( aa_sequences.GetSize()) + " chains.");

        // skip undefined pdbs
        if( aa_sequences.IsEmpty())
        {
          BCL_MessageCrt( "pdb does not contain chains.");
          continue;
        }

        // skip undesired chains, if m_ChainIds is empty, then analyze all chains
        if( !m_ChainIds.empty())
        {
          util::SiPtrVector< const biol::AASequence> aa_sequences_desired;
          for
          (
            util::SiPtrVector< const biol::AASequence>::const_iterator
              seq_itr( aa_sequences.Begin()), seq_itr_end( aa_sequences.End());
            seq_itr != seq_itr_end; ++seq_itr
          )
          {
            if( m_ChainIds.find( ( *seq_itr)->GetChainID()) != std::string::npos)
            {
              aa_sequences_desired.PushBack( *seq_itr);
            }
          }
          aa_sequences.InternalData().swap( aa_sequences_desired.InternalData());
        }

        // iterate through all chains
        for
        (
          util::SiPtrVector< const biol::AASequence>::const_iterator
            seq_itr( aa_sequences.Begin()), seq_itr_end( aa_sequences.End());
          seq_itr != seq_itr_end;
          ++seq_itr
        )
        {
          // iterate over the current chain to get the first amino acid residue
          for
          (
            biol::AASequence::const_iterator
              aa_1_itr( ( *seq_itr)->GetData().Begin()), aa_itr_end( ( *seq_itr)->GetData().End());
            aa_1_itr != aa_itr_end;
            ++aa_1_itr
          )
          {
            if( !( *aa_1_itr)->GetType().IsDefined())
            {
              continue;
            }
            // proceed only if coordinates are given and it is a natural amino acid
            if( !( *aa_1_itr)->GetType()->IsNaturalAminoAcid())
            {
              continue;
            }
            const linal::Vector3D &cb_a( ( *aa_1_itr)->GetFirstSidechainAtom().GetCoordinates());
            const linal::Vector3D &ca_a( ( *aa_1_itr)->GetCA().GetCoordinates());
            if( !ca_a.IsDefined() || !cb_a.IsDefined())
            {
              continue;
            }
            const util::SiPtr< const assemble::SSE> sse1( ( *protein_model_itr)->GetSSE( **aa_1_itr));
            ++number_aas( ( *aa_1_itr)->GetType().GetIndex());

            // iterate over the current chain and all other chains
            for( util::SiPtrVector< const biol::AASequence>::const_iterator seq_2_itr( seq_itr); seq_2_itr != seq_itr_end; ++seq_2_itr)
            {
              // iterate over the current chain to get the second amino acid residue to which the distance of the first amino acid residue is calculated
              for
              (
                biol::AASequence::const_iterator
                  aa_2_itr( seq_2_itr == seq_itr ? aa_1_itr + 1 : ( *seq_2_itr)->GetData().Begin());
                  aa_2_itr < ( *seq_2_itr)->GetData().End();
                ++aa_2_itr
              )
              {
                if( !( *aa_2_itr)->GetType().IsDefined())
                {
                  continue;
                }
                // proceed only if coordinates are given and it is a natural amino acid
                if( !( *aa_2_itr)->GetFirstSidechainAtom().GetCoordinates().IsDefined()
                    || !( *aa_2_itr)->GetType()->IsNaturalAminoAcid())
                {
                  continue;
                }
                const linal::Vector3D &cb_b( ( *aa_2_itr)->GetFirstSidechainAtom().GetCoordinates());
                const linal::Vector3D &ca_b( ( *aa_2_itr)->GetCA().GetCoordinates());
                if( !ca_b.IsDefined() || !cb_b.IsDefined())
                {
                  continue;
                }

                // get sequence separation
                const size_t aa_aa_seq_separation( biol::SequenceSeparation( ( **aa_1_itr), ( **aa_2_itr)));

                // skip amino acid residues that are not from the current chain or are too close
                if( !util::IsDefined( aa_aa_seq_separation) || aa_aa_seq_separation <= m_AADistSeqExcl)
                {
                  continue;
                }

                // skip atoms on the same SSE
                util::SiPtr< const assemble::SSE> sse2( ( *protein_model_itr)->GetSSE( **aa_2_itr));
                if( sse1 == sse2 && sse1.IsDefined())
                {
                  continue;
                }
                ++number_pairs;

                const double cb_cb_distance( linal::Distance( cb_a, cb_b));
                // distance weight
                const double nearest_neighbor_separation
                (
                  cb_cb_distance > 12.0 ? 2.0 : biol::NearestAtomVdWSphereSeparation( **aa_1_itr, **aa_2_itr, true)
                );
                if
                (
                  !nearest_neighbor_separation
                  &&
                  (
                    ( *aa_1_itr)->GetType() != biol::GetAATypes().CYS
                    || ( *aa_2_itr)->GetType() != biol::GetAATypes().CYS
                  )
                )
                {
                  // skip clashes except disulfide bonds
                  BCL_MessageStd
                  (
                    "Clash between " + ( *aa_1_itr)->GetIdentification()
                    + " and " + ( *aa_2_itr)->GetIdentification()
                    + " in " + model_name
                  );
                  continue;
                }

                if
                (
                  m_VdwRadiusStats
                  && cb_cb_distance < 8.0
                  && (
                       ( *aa_1_itr)->GetType() != biol::GetAATypes().CYS
                       || ( *aa_2_itr)->GetType() != biol::GetAATypes().CYS
                     )
                )
                {
                  // iterate over all atom pairs
                  for
                  (
                    auto itr_atom_a( ( *aa_1_itr)->GetAtoms().Begin()), itr_atom_a_end( ( *aa_1_itr)->GetAtoms().End());
                    itr_atom_a != itr_atom_a_end;
                    ++itr_atom_a
                  )
                  {
                    if
                    (
                      ( ( *itr_atom_a)->GetType()->IsBackBone()
                        || ( *itr_atom_a)->GetType() == ( *aa_1_itr)->GetType()->GetFirstSidechainAtomType())
                      && ( *itr_atom_a)->GetCoordinates().IsDefined()
                    )
                    {
                      const std::string atom_name_a
                      (
                        ( *itr_atom_a)->GetType()->IsBackBone()
                        ? ( *itr_atom_a)->GetType()->GetName()
                        : ( *aa_1_itr)->GetType()->GetName()
                      );
                      const double vdw_radii_a( GetVdwRadiiOfType( atom_name_a));
                      for
                      (
                        auto itr_atom_b( ( *aa_2_itr)->GetAtoms().Begin()), itr_atom_b_end( ( *aa_2_itr)->GetAtoms().End());
                        itr_atom_b != itr_atom_b_end;
                        ++itr_atom_b
                      )
                      {
                        if
                        (
                          ( ( *itr_atom_b)->GetType()->IsBackBone()
                            || ( *itr_atom_b)->GetType() == ( *aa_2_itr)->GetType()->GetFirstSidechainAtomType())
                          && ( *itr_atom_b)->GetCoordinates().IsDefined()
                        )
                        {
                          const std::string atom_name_b
                          (
                            ( *itr_atom_b)->GetType()->IsBackBone()
                            ? ( *itr_atom_b)->GetType()->GetName()
                            : ( *aa_2_itr)->GetType()->GetName()
                          );
                          const double vdw_radii_b( GetVdwRadiiOfType( atom_name_b));
                          const double dist( linal::Distance( ( *itr_atom_b)->GetCoordinates(), ( *itr_atom_a)->GetCoordinates()));
                          if( dist >= vdw_radii_a + vdw_radii_b)
                          {
                            if( atom_name_a < atom_name_b)
                            {
                              if( !util::StartsWith( atom_name_a, "GLYCINE") || !util::StartsWith( atom_name_b, "O"))
                              {
                                vdw_map[ atom_name_a + " " + atom_name_b].PushBack( dist);
                              }
                            }
                            else
                            {
                              vdw_map[ atom_name_b + " " + atom_name_a].PushBack( dist);
                            }
                          }
                          else
                          {
                            BCL_MessageStd
                            (
                              "Partial clash: " + ( *aa_1_itr)->GetIdentification() + " " + atom_name_a
                              + ( *aa_2_itr)->GetIdentification() + " " + atom_name_b
                            );
                          }
                        }
                      }
                    }
                  }
                }

                fraction_beyond_12A( ( *aa_1_itr)->GetType().GetIndex())( ( *aa_2_itr)->GetType().GetIndex()) +=
                  cb_cb_distance > 12.0 ? 1.0 : 0.0;
                fraction_beyond_12A( ( *aa_2_itr)->GetType().GetIndex())( ( *aa_1_itr)->GetType().GetIndex()) +=
                  cb_cb_distance > 12.0 ? 1.0 : 0.0;
                const double aa_aa_distance_weight( nearest_neighbor_separation < 1.0 ? 1.0 : 0.0);
                if( aa_aa_distance_weight)
                {
                  ++number_interactions( ( *aa_1_itr)->GetType().GetIndex());
                  ++number_interactions( ( *aa_2_itr)->GetType().GetIndex());
                }
                if( cb_cb_distance > 12.0)
                {
                  continue;
                }
                //
                // add amino acid pair distance to the map
                biol::AAType first_type( std::min( ( *aa_1_itr)->GetType(), ( *aa_2_itr)->GetType()));
                biol::AAType second_type( std::max( ( *aa_1_itr)->GetType(), ( *aa_2_itr)->GetType()));

                double angle( linal::ProjAngleCosinus( linal::UnitVector( ca_a, cb_a), linal::UnitVector( ca_a, cb_b)));
                double angleb( linal::ProjAngleCosinus( linal::UnitVector( ca_b, cb_b), linal::UnitVector( ca_b, cb_a)));
                biol::SSType sstype_1( sse1.IsDefined() ? sse1->GetType() : biol::GetSSTypes().COIL);
                biol::SSType sstype_2( sse2.IsDefined() ? sse2->GetType() : biol::GetSSTypes().COIL);
                if( ( *aa_1_itr)->GetType() > ( *aa_2_itr)->GetType())
                {
                  std::swap( angle, angleb);
                  std::swap( sstype_1, sstype_2);
                }
                std::pair< biol::AAType, biol::AAType> types( first_type, second_type);
                //BCL_MessageStd( "Scorestat: " + first_type.GetName() + " " + second_type.GetName() + " " + util::Format()( cb_cb_distance) + " " + util::Format()( angleb) + " " + util::Format()( angle));
                distance_map_background[ sstype_1][ sstype_2][ types]->PushBack( cb_cb_distance, angleb, angle);
                distance_map_contacts[ types]->PushBack( cb_cb_distance, angleb, angle, aa_aa_distance_weight);
                if( first_type == second_type)
                {
                  distance_map_background[ sstype_2][ sstype_1][ types]->PushBack( cb_cb_distance, angle, angleb);
                  distance_map_contacts[ types]->PushBack( cb_cb_distance, angle, angleb, aa_aa_distance_weight);
                }
              } // end of iterating the current chain to get the second amino acid residues
            } // end of iterating the current chain and all other chains
          } // end of iteration over the current chain to get the first amino acid residues
        } // end of iteration over all chains
      } // end of iteration over all models

      math::RunningAverage< double> ave_frac_beyond_12A;
      for( auto itr_f( fraction_beyond_12A.Begin()), itr_f_end( fraction_beyond_12A.End()); itr_f != itr_f_end; ++itr_f)
      {
        for( auto itr_g( itr_f->Begin()), itr_g_end( itr_f->End()); itr_g != itr_g_end; ++itr_g)
        {
          ave_frac_beyond_12A.AddWeightedObservation( itr_g->GetAverage(), itr_g->GetWeight());
        }
      }
      linal::Matrix< float> energy_beyond_12A
      (
        biol::AATypes::s_NumberStandardAATypes,
        biol::AATypes::s_NumberStandardAATypes
      );
      for( size_t i( 0); i < biol::AATypes::s_NumberStandardAATypes; ++i)
      {
        for( size_t j( 0); j < biol::AATypes::s_NumberStandardAATypes; ++j)
        {
          energy_beyond_12A( i, j) = -log( fraction_beyond_12A( i)( j) / ave_frac_beyond_12A.GetAverage());
        }
      }

      for
      (
        auto itr_ss_type_outer( distance_map_background.begin()), itr_ss_type_outer_end( distance_map_background.end());
        itr_ss_type_outer != itr_ss_type_outer_end;
        ++itr_ss_type_outer
      )
      {
        for
        (
          auto itr_ss_type_inner( itr_ss_type_outer->begin()), itr_ss_type_inner_end( itr_ss_type_outer->end());
          itr_ss_type_inner != itr_ss_type_inner_end;
          ++itr_ss_type_inner
        )
        {
          for
          (
            auto itr_bg( itr_ss_type_inner->begin()),
                 itr_bg_end( itr_ss_type_inner->end()),
                 itr_sim_bg( rand_bg.begin());
            itr_bg != itr_bg_end;
            ++itr_bg, ++itr_sim_bg
          )
          {
            math::Tensor< double> &bg_hist( itr_bg->second->GetChangeableHistogram());
            const math::Tensor< double> &sim_hist( itr_sim_bg->second->GetHistogram());
            const double total_pseudocount_bg( 32.0 * bg_hist.GetValues().GetSize());

            const double sim_hist_raw_sum( sim_hist.GetValues().Sum());
            const double bg_hist_raw_sum( bg_hist.GetValues().Sum());

            const double bg_hist_sum( bg_hist_raw_sum + total_pseudocount_bg);
            const double total_pseudocount_sim( total_pseudocount_bg * ( sim_hist_raw_sum / bg_hist_raw_sum));
            const double sim_hist_sum( sim_hist_raw_sum + total_pseudocount_sim);
            const double pseudocount_sim( total_pseudocount_sim / double( bg_hist.GetValues().GetSize()));
            const double size_ratio( sim_hist_sum / bg_hist_sum);
            //        BCL_MessageStd
            //        (
            //          "PC: " + itr_bg->first.first.GetName()
            //          + " " + itr_bg->first.second.GetName() + " "
            //          + util::Format()( bg_hist_raw_sum) + " "
            //          + util::Format()( sim_hist_raw_sum) + " "
            //          + util::Format()( pseudocount_sim) + " "
            //          + util::Format()( size_ratio) + " "
            //          + util::Format()( bg_hist_sum) + " "
            //          + util::Format()( sim_hist_sum)
            //        );
            auto itr_sim_bge( sim_hist.Begin());
            const double twelve_a_propensity( -energy_beyond_12A( itr_bg->first.first, itr_bg->first.second));
            for
            (
              auto itr_bge( bg_hist.Begin()), itr_bge_end( bg_hist.End());
              itr_bge != itr_bge_end;
              ++itr_bge, ++itr_sim_bge
            )
            {
              *itr_bge = -log( ( *itr_bge + 32.0) * size_ratio / ( *itr_sim_bge + pseudocount_sim)) + twelve_a_propensity;
            }
          }
        }
      }

      storage::Map< std::string, double> min_vdw_radii;
      storage::Map< std::string, double> vdw_radii_final;
      if( m_VdwRadiusStats)
      {
        for( auto itr( vdw_map.Begin()), itr_end( vdw_map.End()); itr != itr_end; ++itr)
        {
          auto strings( util::SplitString( itr->first, " "));
          const std::string &str_a( strings( 0)), &str_b( strings( 1));
          min_vdw_radii[ str_a] = GetVdwRadiiOfType( str_a);
          min_vdw_radii[ str_b] = GetVdwRadiiOfType( str_b);
        }
        vdw_radii_final = chemistry::BondLengths::ComputeVdwRadii( vdw_map, min_vdw_radii);
      }

      // write statistics
      std::ostringstream stream;

      for( size_t i( 0), n_ss_types( 3); i < n_ss_types; ++i)
      {
        for( size_t j( 0), n_ss_types( 3); j < n_ss_types; ++j)
        {
          // get the first residue
          for
          (
            biol::AATypes::const_iterator aa_1_itr( biol::GetAATypes().Begin()), aa_itr_end( biol::GetAATypes().End());
              aa_1_itr != aa_itr_end;
            ++aa_1_itr
          )
          {
            // get the second residue
            for
            (
              biol::AATypes::const_iterator aa_2_itr( aa_1_itr);
                aa_2_itr != aa_itr_end;
              ++aa_2_itr
            )
            {
              const std::pair< biol::AAType, biol::AAType> aa_pair( *aa_1_itr, *aa_2_itr);

              // check to see if distance_map_contacts has the entry
              const std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram3D> >::const_iterator
                distance_map_contacts_entry( distance_map_background[ i][ j].find( aa_pair));

              if( distance_map_contacts_entry != distance_map_background[ i][ j].end())
              {
                stream << ( *aa_1_itr)->GetOneLetterCode() << " " << ( *aa_2_itr)->GetOneLetterCode() << "\n";
                stream << *distance_map_contacts_entry->second << '\n';
              }
            }
          }
        }
      }

      BCL_MessageStd
      (
        "# Pairs: " + util::Format()( number_pairs)
        + "\ntotal interactions: " + util::Format()( number_interactions.Sum())
        + "\n# aas: " + util::Format()( number_aas.Sum())
      );
      BCL_MessageStd
      (
        "# P(Pair interacting): " + util::Format()( double( number_interactions.Sum()) / double( number_pairs))
        + "\nP(AA interacting with any): " + util::Format()( double( number_interactions.Sum()) / double( number_aas.Sum()))
      );
      BCL_MessageStd( "AAType\t# of aa type\t# of interactions w AAType\tFraction of AATypes interacting");
      for( size_t i( 0); i < size_t( 20); ++i)
      {
        biol::AAType type( i);
        BCL_MessageStd( std::string( 1, type->GetOneLetterCode()) + "\t" + util::Format()( number_aas( i)) + "\t" + util::Format()( number_interactions( i)) + "\t" + util::Format()( double( number_interactions( i)) / double( number_aas( i))));
      }

      return stream.str();

    } // end of operator

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AADistanceAngleContacts::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Computes amino acid pair distance statistics.");

      parameters.AddInitializer
      (
        "min count",
        "Minimum count per bin to not consider it a clashing bin",
        io::Serialization::GetAgentWithMin( &m_MinCounts, 0.0),
        "10"
      );

      parameters.AddInitializer
      (
        "bin_size",
        "the bin size for the histogram",
        io::Serialization::GetAgent( &m_BinSize),
        "1"
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "seq_exlusion",
        "sequence exclusion: sequence separation below which aa distance is not considered",
        io::Serialization::GetAgent( &m_AADistSeqExcl),
        "2"
      );
      parameters.AddInitializer
      (
        "vdw_stats",
        "Whether to collect and VdW radius stats for individual atom types. "
        "This modestly increases overall time and memory consumption by a few percent",
        io::Serialization::GetAgent( &m_VdwRadiusStats),
        "false"
      );
      return parameters;
    }
  } // namespace scorestat
} // namespace bcl

