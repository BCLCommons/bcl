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

// include header for this class
#include "scorestat/bcl_scorestat_aa_count.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AACount::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AACount())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &AACount::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Table",
        GetStaticClassName< AACount::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &AACount::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "aa_count.tbl",
        GetStaticClassName< AACount::OutputOption>()
      };
      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AACount::AACount():
      m_OutputOption( e_Table),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    AACount *AACount::Clone() const
    {
      return new AACount( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &AACount::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AACount::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &AACount::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AACount::GetAlias() const
    {
      static std::string s_name( "AACount");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AACount::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // aa count distribution according to environmental types, such as membrane core, transition, soluble
      storage::Vector< storage::Vector< size_t> > aa_region_count
      (
        biol::GetAATypes().GetEnumCount(),
        storage::Vector< size_t>( biol::GetEnvironmentTypes().GetEnumCount(), size_t( 0))
      );

      // aa count distribution according to secondary structure elements and environmental types
      storage::Vector< storage::Vector< storage::Vector< size_t> > > aa_sse_region_count
      (
        biol::GetSSTypes().COIL.GetIndex() + 1,
        storage::Vector< storage::Vector< size_t> >
        (
          biol::GetAATypes().GetEnumCount(),
          storage::Vector< size_t>( biol::GetEnvironmentTypes().GetEnumCount(), size_t( 0))
        )
      );

      // iterate through all protein models
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get the current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get protein pdb filename
        const util::ShPtr< util::Wrapper< std::string> >
          &sp_model_filename( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_filename( sp_model_filename.IsDefined() ? sp_model_filename->GetData() : "");

        // get all chains in current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // get membrane for current protein model
        const util::ShPtr< biol::Membrane> sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        if( !sp_membrane.IsDefined())
        {
          BCL_MessageDbg( util::Format()( model_filename) + " does not have an associated membrane, assuming it's a soluble protein.");
        }

        // iterate over all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator
            chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // get all sses in current chain
          const storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>
            &all_sses( ( *chain_itr)->GetData());

          // iterate over all sses in current chain
          for
          (
            storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
              sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            // skip protein models without sse entries in the pdb file
            if( ( *chain_itr)->GetNumberSSEs() == 1 && ( *sse_itr)->GetType() == biol::GetSSTypes().COIL)
            {
              BCL_MessageStd( util::Format()( model_filename) + " probably does not have sse entries in the pdb file, skipping");
              continue;
            }

            // skip undefined sses
            if( !( *sse_itr)->GetType().IsDefined())
            {
              continue;
            }

            // get the type of current sse
            const biol::SSType &current_sse_type( ( *sse_itr)->GetType());

            // get all amino acids in current sse
            const util::ShPtrVector< biol::AABase> &all_amino_acids( ( *sse_itr)->GetData());

            // iterate over all amino acids in current sse
            for
            (
              util::ShPtrVector< biol::AABase>::const_iterator
                aa_itr( all_amino_acids.Begin()), aa_itr_end( all_amino_acids.End());
              aa_itr != aa_itr_end;
              ++aa_itr
            )
            {
              // get environment type in which the current amino acid can be found
              const biol::EnvironmentType current_environment_type
              (
                sp_membrane.IsDefined() ?
                    sp_membrane->DetermineEnvironmentType( ( *aa_itr)->GetFirstSidechainAtom().GetCoordinates()) :
                      biol::GetEnvironmentTypes().e_Solution
              );

              // skip undefined environment
              if( !current_environment_type.IsDefined())
              {
                continue;
              }

              ++aa_region_count( ( *aa_itr)->GetType())( current_environment_type);
              ++aa_sse_region_count( current_sse_type)( ( *aa_itr)->GetType())( current_environment_type);
            } // end of iterating over all amino acids
          } // end of iterating over all sses in current chain
        } // end of iterating over all chains in current protein model
      } // end of iterating over all protein models in the ensemble

      // write statistics
      std::stringstream stream;

      if( m_OutputOption == e_Table)
      {
        stream << "# aa count in each environment type" << '\n';
        stream << '\n';
        stream << "environment_type" << '\t';

        for
        (
          biol::AATypes::const_iterator aa_itr( biol::GetAATypes().Begin()),
            aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          stream << ( *aa_itr)->GetOneLetterCode() << '\t';
        }
        stream << '\n';

        for
        (
          biol::EnvironmentTypes::const_iterator env_itr( biol::GetEnvironmentTypes().Begin()),
            env_itr_end( biol::GetEnvironmentTypes().End());
          env_itr != env_itr_end;
          ++env_itr
        )
        {
          stream << *env_itr << '\t';

          // all amino acid types
          for
          (
            biol::AATypes::const_iterator aa_itr( biol::GetAATypes().Begin()),
              aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
            aa_itr != aa_itr_end;
            ++aa_itr
          )
          {
            stream << aa_region_count( *aa_itr)( *env_itr) << '\t';
          }
          stream << '\n';
        }

        stream << '\n';
        stream << "# aa count in each environment type categorized by sse type" << '\n';

        for
        (
          biol::SSTypes::const_iterator sse_itr( biol::GetSSTypes().Begin()),
            sse_itr_end( biol::GetSSTypes().COIL.GetIterator() + 1);
          sse_itr != sse_itr_end;
          ++sse_itr
        )
        {
          stream << '\n';
          stream << ( *sse_itr)->GetName() << '\n';
          stream << "environment_type" << '\t';

          for
          (
            biol::AATypes::const_iterator aa_itr( biol::GetAATypes().Begin()),
              aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
            aa_itr != aa_itr_end;
            ++aa_itr
          )
          {
            stream << ( *aa_itr)->GetOneLetterCode() << '\t';
          }
          stream << '\n';

          for
          (
            biol::EnvironmentTypes::const_iterator env_itr( biol::GetEnvironmentTypes().Begin()),
              env_itr_end( biol::GetEnvironmentTypes().End());
            env_itr != env_itr_end;
            ++env_itr
          )
          {
            stream << *env_itr << '\t';

            // all amino acid types
            for
            (
              biol::AATypes::const_iterator aa_itr( biol::GetAATypes().Begin()),
                aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
              aa_itr != aa_itr_end;
              ++aa_itr
            )
            {
              stream << aa_sse_region_count( *sse_itr)( *aa_itr)( *env_itr) << '\t';
            }
            stream << '\n';
          }
        }
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AACount::GetSerializer() const
    {

      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes amino acid occurrence statistics."
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "Table"
      );
      return parameters;
    }
  } // namespace scorestat
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
#include "scorestat/bcl_scorestat_aa_distance.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score_aa_pair_sidechain_interaction.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AADistance::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AADistance())
    );

    //! @brief CategoryOptionName as string
    //! @param CATEGORY_OPTION the AADistanceCategoryOption
    //! @return the string for the CategoryOptionName
    const std::string &AADistance::GetCategoryOptionName( const AADistanceCategoryOption &CATEGORY_OPTION)
    {
      static const std::string s_names[] =
      {
          "OneChainNew",
          "OneChainOld",
          "AllChainNew",
          "AllChainOld",
          GetStaticClassName< AADistanceCategoryOption>()
      };
      return s_names[ CATEGORY_OPTION];
    }

    //! @brief CategoryFileName as string
    //! @param CATEGORY_OPTION the AADistanceCategoryOption
    //! @return the string for the CategoryFileName
    const std::string &AADistance::GetCategoryFileName( const AADistanceCategoryOption &CATEGORY_OPTION)
    {
      static const std::string s_names[] =
      {
          "aa_distances_one_chain_new.histograms",
          "aa_distances_one_chain_old.histograms",
          "aa_distances_all_chain_new.histograms",
          "aa_distances_all_chain_old.histograms",
          GetStaticClassName< AADistanceCategoryOption>()
      };
      return s_names[ CATEGORY_OPTION];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AADistance::AADistance() :
      m_CategoryOption( e_OneChainNew),
      m_BinSize( 1.0),
      m_ChainIds( ""),
      m_AADistSeqExcl( 6)
    {
    }

    //! @brief virtual copy constructor
    AADistance *AADistance::Clone() const
    {
      return new AADistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AADistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AADistance::GetOutFilePostfix() const
    {
      return GetCategoryFileName( m_CategoryOption);
    }

    //! @brief returns the binsize for the histogram
    //! @return the binsize for the histogram
    const double &AADistance::GetBinSize() const
    {
      return m_BinSize;
    }

    //! @brief returns chain id
    //! @return chain id
    const std::string &AADistance::GetChainId() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AADistance::GetAlias() const
    {
      static const std::string s_name( "AADistance");
      return s_name;
    }

  //////////////
  // operator //
  //////////////

    //! @brief add statistics about the given protein and chains to the internal table and histograms
    //! @param ENSEMBLE the protein ensemble the statistics of which to be computed
    std::string AADistance::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // initializes maps for storing histograms of distances between all types of amino acid pairs
      std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> > distance_map_one_old;
      std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> > distance_map_one_new;
      std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> > distance_map_all_old;
      std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> > distance_map_all_new;

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

          // initialize a histogram
          util::ShPtr< math::Histogram> histogram(
            new math::Histogram( double( 0), m_BinSize, size_t( 100 / m_BinSize)));

          // add aa pair-histogram to the map as a side effect of subscripting the map
          distance_map_one_old[ std::pair< biol::AAType, biol::AAType>( *aa_1_itr, *aa_2_itr)] = histogram.HardCopy();
          distance_map_one_new[ std::pair< biol::AAType, biol::AAType>( *aa_1_itr, *aa_2_itr)] = histogram.HardCopy();
          distance_map_all_old[ std::pair< biol::AAType, biol::AAType>( *aa_1_itr, *aa_2_itr)] = histogram.HardCopy();
          distance_map_all_new[ std::pair< biol::AAType, biol::AAType>( *aa_1_itr, *aa_2_itr)] = histogram.HardCopy();
        }
      }

      // initialize histogram for storing disulfide bond statistics
      storage::VectorND< 2, math::Histogram> disulfide_bonds
      ( math::Histogram( double( 0), m_BinSize, size_t( 100 / m_BinSize)));

      // cutoff of disulfide bond length
      static const double disulfide_bond_cutoff( 2.5);

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
        BCL_MessageStd( "building sequences from pdb chains");
        util::SiPtrVector< const biol::AASequence> aa_sequences( protein_model.GetSequences());

        // number of chains
        BCL_MessageStd( "pdb has " + util::Format()( aa_sequences.GetSize()) + " chains.");

        // skip undefined pdbs
        if( aa_sequences.IsEmpty())
        {
          BCL_MessageCrt( "pdb does not contain chains.");
          continue;
        }

        // iterate through all chains
        for
        (
          util::SiPtrVector< const biol::AASequence>::const_iterator
            seq_itr( aa_sequences.Begin()), seq_itr_end( aa_sequences.End());
          seq_itr != seq_itr_end; ++seq_itr
        )
        {
          // skip undesired chains, if m_ChainIds is empty, then analyze all chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *seq_itr)->GetChainID()) == std::string::npos)
          {
            BCL_MessageDbg( std::string( "Skip chain: ") + ( *seq_itr)->GetChainID() + std::string( ", chain ids to use: ") + m_ChainIds);

            // continue to the next chain
            continue;
          }

          // iterate over the current chain to get the first amino acid residue
          for
          (
            biol::AASequence::const_iterator
              aa_1_itr( ( *seq_itr)->GetData().Begin()), aa_itr_end( ( *seq_itr)->GetData().End());
            aa_1_itr != aa_itr_end;
            ++aa_1_itr
          )
          {
            // proceed only if coordinates are given and it is a natural amino acid
            if
            (
              !( *aa_1_itr)->GetFirstSidechainAtom().GetCoordinates().IsDefined() ||
              !( *aa_1_itr)->GetType().IsDefined() || !( *aa_1_itr)->GetType()->IsNaturalAminoAcid()
            )
            {
              continue;
            }

            // iterate over the current chain and all other chains
            for( util::SiPtrVector< const biol::AASequence>::const_iterator seq_2_itr( seq_itr); seq_2_itr != seq_itr_end; ++seq_2_itr)
            {
              // iterate over the current chain to get the second amino acid residue to which the distance of the first amino acid residue is calculated
              for
              (
                biol::AASequence::const_iterator aa_2_itr( ( *seq_itr)->GetData().Begin() + 1);
                  aa_2_itr != ( *seq_itr)->GetData().End();
                ++aa_2_itr
              )
              {
                // proceed only if coordinates are given and it is a natural amino acid
                if( !( *aa_2_itr)->GetFirstSidechainAtom().GetCoordinates().IsDefined()
                    || !( *aa_2_itr)->GetType().IsDefined() || !( *aa_2_itr)->GetType()->IsNaturalAminoAcid())
                {
                  continue;
                }

                // check for disulfide bond
                if( ( *aa_1_itr)->GetType() == biol::GetAATypes().CYS && ( *aa_2_itr)->GetType() == biol::GetAATypes().CYS)
                {
                  // calculate distance between two cysteine residues
                  const double cys_cys_distance
                  (
                    biol::Distance( ( *aa_1_itr)->GetAtom( biol::GetAtomTypes().SG), ( *aa_2_itr)->GetAtom( biol::GetAtomTypes().SG))
                  );

                  // store the cys_cys_distance if it is defined and is within the range of disulfide bond length
                  if( util::IsDefined( cys_cys_distance) && cys_cys_distance < disulfide_bond_cutoff)
                  {
                    // aa pair distance
                    const double aa_aa_distance( biol::FirstSidechainAtomDistance( ( **aa_1_itr), ( **aa_2_itr)));

                    // insert into histogram
                    disulfide_bonds( 0).PushBack( aa_aa_distance);
                    disulfide_bonds( 1).PushBack( cys_cys_distance);

                    // skip the current amino acid residue since its distance has been calculated
                    continue;
                  }
                }

                // get sequence separation
                const size_t aa_aa_seq_separation( biol::SequenceSeparation( ( **aa_1_itr), ( **aa_2_itr)));

                // skip amino acid residues that are not from the current chain or are too close
                if( !util::IsDefined( aa_aa_seq_separation) || aa_aa_seq_separation < m_AADistSeqExcl)
                {
                  continue;
                }

                // distance between the first side chain atom of the first amino acid and that of the second amino acid
                const double aa_aa_distance( biol::FirstSidechainAtomDistance( ( **aa_1_itr), ( **aa_2_itr)));

                // distance weight
                const double aa_aa_distance_weight
                (
                  score::AAPairSidechainInteraction::WeightOfInteraction( ( **aa_1_itr), ( **aa_2_itr))
                );

                // add amino acid pair distance to the map
                biol::AAType first_type( std::min( ( *aa_1_itr)->GetType(), ( *aa_2_itr)->GetType()));
                biol::AAType second_type( std::max( ( *aa_1_itr)->GetType(), ( *aa_2_itr)->GetType()));
                if( seq_itr == seq_2_itr)
                {
                  distance_map_one_old[ std::pair< biol::AAType, biol::AAType>( first_type, second_type)]->PushBack( aa_aa_distance);
                  distance_map_one_new[ std::pair< biol::AAType, biol::AAType>( first_type, second_type)]->PushBack( aa_aa_distance, aa_aa_distance_weight);
                }
                else
                {
                  distance_map_all_old[ std::pair< biol::AAType, biol::AAType>( first_type, second_type)]->PushBack( aa_aa_distance);
                  distance_map_all_new[ std::pair< biol::AAType, biol::AAType>( first_type, second_type)]->PushBack( aa_aa_distance, aa_aa_distance_weight);
                }
              } // end of iterating the current chain to get the second amino acid residues
            } // end of iterating the current chain and all other chains
          } // end of iteration over the current chain to get the first amino acid residues
        } // end of iteration over all chains
      } // end of iteration over all models

      // write statistics
      std::ostringstream stream;

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

          if( m_CategoryOption == e_AllChainNew)
          {
            // check to see if distance_map_all_new has the entry
            const std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> >::const_iterator
              distance_map_all_chain_new_entry( distance_map_all_new.find( aa_pair));

            if( distance_map_all_chain_new_entry != distance_map_all_new.end() && !distance_map_all_chain_new_entry->second->IsEmpty())
            {
              stream << ( *aa_1_itr)->GetOneLetterCode() << " " << ( *aa_2_itr)->GetOneLetterCode() << "\n";
              stream << *distance_map_all_chain_new_entry->second;
            }
          }
          else if( m_CategoryOption == e_AllChainOld)
          {
            // check to see if distance_map_all_old has the entry
            const std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> >::const_iterator
              distance_map_all_old_entry( distance_map_all_old.find( aa_pair));

            if( distance_map_all_old_entry != distance_map_all_old.end() && !distance_map_all_old_entry->second->IsEmpty())
            {
              stream << ( *aa_1_itr)->GetOneLetterCode() << " " << ( *aa_2_itr)->GetOneLetterCode() << "\n";
              stream << *distance_map_all_old_entry->second;
            }
          }
          else if( m_CategoryOption == e_OneChainNew)
          {
            // check to see if distance_map_one_new has the entry
            const std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> >::const_iterator
              distance_map_one_new_entry( distance_map_one_new.find( aa_pair));

            if( distance_map_one_new_entry != distance_map_one_new.end() && !distance_map_one_new_entry->second->IsEmpty())
            {
              stream << ( *aa_1_itr)->GetOneLetterCode() << " " << ( *aa_2_itr)->GetOneLetterCode() << "\n";
              stream << *distance_map_one_new_entry->second;
            }
          }
          else if( m_CategoryOption == e_OneChainOld)
          {
            // check to see if distance_map_one_old has the entry
            const std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> >::const_iterator
              distance_map_one_old_entry( distance_map_one_old.find( aa_pair));

            if( distance_map_one_old_entry != distance_map_one_old.end() && !distance_map_one_old_entry->second->IsEmpty())
            {
              stream << ( *aa_1_itr)->GetOneLetterCode() << " " << ( *aa_2_itr)->GetOneLetterCode() << "\n";
              stream << *distance_map_one_old_entry->second;
            }
          }
        }
      }

      return stream.str();

    } // end of operator

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AADistance::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Computes amino acid pair distance statistics.");

      parameters.AddInitializer
      (
        "category",
        "the distance category: OneChainNew, OneChainOld, AllChainNew, AllChainOld",
        io::Serialization::GetAgent( &m_CategoryOption),
        "OneChainNew"
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
        "6"
      );

      return parameters;
    }
  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_aa_distance_matrix.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_table.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    //! single instance of this class
    const util::SiPtr< util::ObjectInterface> AADistanceMatrix::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AADistanceMatrix())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &AADistanceMatrix::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_names[] =
      {
          "Table",
          GetStaticClassName< AADistanceMatrix::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &AADistanceMatrix::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_names[] =
      {
        "aa_distance_matrix.tbl",
        GetStaticClassName< AADistanceMatrix::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

      //! @brief implicit constructor
      AADistanceMatrix::AADistanceMatrix() :
          m_OutputOption( e_Table),
          m_ColumnChainId( ""),
          m_RowChainId( "")
      {
        // nothing else to do
      }

      //! @brief virtual copy constructor
      //! @return a pointer to an instance of the copied AADistanceMatrix
      AADistanceMatrix *AADistanceMatrix::Clone() const
      {
        return new AADistanceMatrix( *this);
      }

  /////////////////
  // data access //
  /////////////////

      //! @brief returns class name
      //! @return the class name as constant reference to std::string
      const std::string &AADistanceMatrix::GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return the string to append to the the end of a filename to identify this analysis
      const std::string &AADistanceMatrix::GetOutFilePostfix() const
      {
        return GetOutputFileName( m_OutputOption);
      }

      //! @brief returns chain ids
      //! @return chain ids
      const std::string &AADistanceMatrix::GetColumnChainId() const
      {
        return m_ColumnChainId;
      }

      //! @brief returns chain ids
      //! @return chain ids
      const std::string &AADistanceMatrix::GetRowChainId() const
      {
        return m_RowChainId;
      }

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &AADistanceMatrix::GetAlias() const
      {
        static const std::string s_name( "AADistanceMatrix");
        return s_name;
      }

  ///////////////
  // operators //
  ///////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the protein ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string AADistanceMatrix::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
      {
        // make sure the protein ensemble is not empty
        BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

        // table header of the matrix for storing aapair distances
        storage::Table< double> distance_matrix;
        storage::Vector< std::string> matrix_column_names;
        std::map< std::pair< std::string, std::string>, double> distance_map;

        // get first protein model from the ensemble
        const assemble::ProteinModel first_protein_model( ( **ENSEMBLE.Begin()));

        // get the desired row chain
        const util::ShPtr< biol::AASequence> row_chain( first_protein_model.GetChain( m_RowChainId[ 0])->GetSequence());

        // get the desired column chain
        const util::ShPtr< biol::AASequence> column_chain( first_protein_model.GetChain( m_ColumnChainId[ 0])->GetSequence());

        // make sure protein model is qualified
        BCL_Assert
        (
          row_chain.IsDefined() && column_chain.IsDefined(),
          "Unqualified protein model: the current protein model does not have desired chains"
        );

        // iterate through the column chain
        for
        (
          util::ShPtrVector< biol::AABase>::const_iterator
            col_aa_itr( column_chain->Begin()), col_aa_itr_end( column_chain->End());
          col_aa_itr != col_aa_itr_end;
          ++col_aa_itr
        )
        {
          // skip undefined residues
          if( !( *col_aa_itr)->GetType()->IsNaturalAminoAcid() || ( *col_aa_itr)->GetAtomCoordinates().IsEmpty())
          {
            continue;
          }
          matrix_column_names.PushBack
          (
            ( *col_aa_itr)->GetType()->GetOneLetterCode() + util::Format()( ( *col_aa_itr)->GetSeqID())
          );
        }

        // initialize the table header
        const util::ShPtr< storage::TableHeader> sp_table_header( new storage::TableHeader( matrix_column_names));
        distance_matrix = storage::Table< double>( sp_table_header);

        // iterate through all protein models in the ensemble
        size_t number_valide_protein_models( 0);
        for
        (
          util::ShPtrVector< assemble::ProteinModel>::const_iterator
            protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
          protein_model_itr != protein_model_itr_end;
          ++protein_model_itr
        )
        {
          // get current protein model
          const assemble::ProteinModel protein_model( **protein_model_itr);

          // get the desired row chain
          const util::ShPtr< biol::AASequence> row_chain( protein_model.GetChain( m_RowChainId[ 0])->GetSequence());

          // get the desired column chain
          const util::ShPtr< biol::AASequence> column_chain( protein_model.GetChain( m_ColumnChainId[ 0])->GetSequence());

          // skip unqualified protein models
          if( !row_chain.IsDefined() || !column_chain.IsDefined())
          {
            BCL_MessageStd( "Unqualified protein model: the current protein model does not have desired chains");
            continue;
          }

          // increment the number of valide protein models
          ++number_valide_protein_models;

          // iterate through row residues to insert rows
          for
          (
            util::ShPtrVector< biol::AABase>::const_iterator
              row_aa_itr( row_chain->Begin()), row_aa_itr_end( row_chain->End());
            row_aa_itr != row_aa_itr_end;
            ++row_aa_itr
          )
          {
            // skip undefined residues
            if( !( *row_aa_itr)->GetType()->IsNaturalAminoAcid() || ( *row_aa_itr)->GetAtomCoordinates().IsEmpty())
            {
              continue;
            }

            // row_aa identifier
            const std::string row_aa_id
            (
              ( *row_aa_itr)->GetType()->GetOneLetterCode() + util::Format()( ( *row_aa_itr)->GetSeqID())
            );

            // iterate through column residues
            for
            (
              util::ShPtrVector< biol::AABase>::const_iterator
                col_aa_itr( column_chain->Begin()), col_aa_itr_end( column_chain->End());
              col_aa_itr != col_aa_itr_end;
              ++col_aa_itr
            )
            {
              // skip undefined residues
              if( !( *col_aa_itr)->GetType()->IsNaturalAminoAcid() || ( *col_aa_itr)->GetAtomCoordinates().IsEmpty())
              {
                continue;
              }

              // col_aa identifier
              const std::string col_aa_id
              (
                ( *col_aa_itr)->GetType()->GetOneLetterCode() + util::Format()( ( *col_aa_itr)->GetSeqID())
              );

              distance_map[ std::pair< std::string, std::string>( row_aa_id, col_aa_id)] +=
                  biol::FirstSidechainAtomDistance( ( **row_aa_itr), ( **col_aa_itr));
            }
          }
        }

        // iterate through row residues to insert rows
        for
        (
          util::ShPtrVector< biol::AABase>::const_iterator
            row_aa_itr( row_chain->Begin()), row_aa_itr_end( row_chain->End());
          row_aa_itr != row_aa_itr_end;
          ++row_aa_itr
        )
        {
          // skip undefined residues
          if( !( *row_aa_itr)->GetType()->IsNaturalAminoAcid() || ( *row_aa_itr)->GetAtomCoordinates().IsEmpty())
          {
            continue;
          }

          // row_aa identifier
          const std::string row_aa_id
          (
            ( *row_aa_itr)->GetType()->GetOneLetterCode() + util::Format()( ( *row_aa_itr)->GetSeqID())
          );

          // create vector for holding distances
          storage::Vector< double> distances;

          // iterate through column residues
          for
          (
            util::ShPtrVector< biol::AABase>::const_iterator
              col_aa_itr( column_chain->Begin()), col_aa_itr_end( column_chain->End());
            col_aa_itr != col_aa_itr_end;
            ++col_aa_itr
          )
          {
            // col_aa identifier
            const std::string col_aa_id
            (
              ( *col_aa_itr)->GetType()->GetOneLetterCode() + util::Format()( ( *col_aa_itr)->GetSeqID())
            );

            distances.PushBack
            (
              distance_map[ std::pair< std::string, std::string>( row_aa_id, col_aa_id)] / number_valide_protein_models
            );
          }

          // insert row to table
          distance_matrix.InsertRow
          (
            row_aa_id,
            distances,
            true
          );

        }
        // write distance matrix
        std::ostringstream stream;
        if( m_OutputOption == e_Table)
        {
          distance_matrix.WriteFormatted( stream);
        }

        return stream.str();
      }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AADistanceMatrix::GetSerializer() const
    {
      io::Serializer paramters;

      paramters.SetClassDescription
      (
        "Compute amino acid distance matrix averaged over a protein conformational ensemble."
      );

      paramters.AddInitializer
      (
        "row_chain_id",
        "id of the chain containing residues corresponding to row names",
        io::Serialization::GetAgent( &m_RowChainId),
        "A"
      );

      paramters.AddInitializer
      (
        "column_chain_id",
        "id of the chain containing residues corresponding to column names",
        io::Serialization::GetAgent( &m_ColumnChainId),
        "A"
      );

      paramters.AddInitializer
      (
        "output",
        "output format",
        io::Serialization::GetAgent( &m_OutputOption),
        "Table"
      );

      return paramters;
    }
  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_contact_order.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score_contact_order.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ContactOrder::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new ContactOrder())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &ContactOrder::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Table",
        "Histogram",
        GetStaticClassName< ContactOrder::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &ContactOrder::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "contact_order.tbl",
        "contact_order.histograms",
        GetStaticClassName< ContactOrder::OutputOption>()
      };
      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

    //! @brief default constructor
    ContactOrder::ContactOrder() :
      m_OutputOption( e_Table),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    ContactOrder *ContactOrder::Clone() const
    {
      return new ContactOrder( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &ContactOrder::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &ContactOrder::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &ContactOrder::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ContactOrder::GetAlias() const
    {
      static std::string s_name( "ContactOrder");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string ContactOrder::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // table for holding contact order statistics
      storage::Table< double> contact_order_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            "nr_aas", "nr_aas_sses", "co_chain_sse", "co_chain_seq", "co_seq", "co_chain_sse_sqr", "co_chain_seq_sqr", "co_seq_sqr"
          )
        )
      );

      // create all types of contact orders
      const contact::Order co_relative_sse( contact::Order::e_RelativeAAsUsed, "co_relative_sse", false);
      const contact::Order co_relative_seq( contact::Order::e_RelativeSequenceLength, "co_relative_length", false);
      const contact::Order co_relative_sse_sqr( contact::Order::e_RelativeSqrAAsUsed, "co_relative_sse_sqr", false);
      const contact::Order co_relative_seq_sqr( contact::Order::e_RelativeSqrSequenceLength, "co_relative_seq_sqr", false);

      // iterate over protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get current protein model name
        const util::ShPtr< util::Wrapper< std::string> >
          &sp_model_name( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( sp_model_name.IsDefined() ? sp_model_name->GetData() : "");

        // get all chains
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // get all sequences
        const util::SiPtrVector< const biol::AASequence> &all_sequences( current_protein_model.GetSequences());

        // iterate over all chains of current protein model
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator
            chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains, if m_ChainIds is empty, analyze all chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // insert rows
          if( !contact_order_table.HasRow( model_name + ( *chain_itr)->GetChainID()))
          {
            contact_order_table.InsertRow( model_name + ( *chain_itr)->GetChainID());
          }

          contact_order_table[ model_name + ( *chain_itr)->GetChainID()][ "nr_aas_sses"] = ( *chain_itr)->GetNumberAAs();
          contact_order_table[ model_name + ( *chain_itr)->GetChainID()][ "co_chain_sse"] = co_relative_sse.ContactOrder( **chain_itr);
          contact_order_table[ model_name + ( *chain_itr)->GetChainID()][ "co_chain_seq"] = co_relative_seq.ContactOrder( **chain_itr);
          contact_order_table[ model_name + ( *chain_itr)->GetChainID()][ "co_chain_sse_sqr"] = co_relative_sse_sqr.ContactOrder( **chain_itr);
          contact_order_table[ model_name + ( *chain_itr)->GetChainID()][ "co_chain_seq_sqr"] = co_relative_seq_sqr.ContactOrder( **chain_itr);
        } // end of iteration over all chains

        // iterate over all sequences
        for
        (
          util::SiPtrVector< const biol::AASequence>::const_iterator
            aaseq_itr( all_sequences.Begin()), aaseq_itr_end( all_sequences.End());
          aaseq_itr != aaseq_itr_end;
          ++aaseq_itr
        )
        {
          if( !contact_order_table.HasRow( model_name + ( *aaseq_itr)->GetChainID()))
          {
            contact_order_table.InsertRow( model_name + ( *aaseq_itr)->GetChainID());
          }

          contact_order_table[ model_name + ( *aaseq_itr)->GetChainID()][ "nr_aas"] = ( *aaseq_itr)->GetSize();
          contact_order_table[ model_name + ( *aaseq_itr)->GetChainID()][ "co_seq"] = co_relative_seq.ContactOrder( **aaseq_itr);
          contact_order_table[ model_name + ( *aaseq_itr)->GetChainID()][ "co_seq_sqr"] = co_relative_seq_sqr.ContactOrder( **aaseq_itr);
        } // end of iteration over sequences
      } // end of iteration over protein ensemble

      // write statistics
      std::stringstream stream;

      if( m_OutputOption == e_Table)
      {
        // contact order summary
        stream << "contact order summary" << '\n';
        contact_order_table.WriteFormatted( stream);
        stream << '\n';
      }
      else if( m_OutputOption == e_Histogram)
      {
        // contact order histograms
        storage::Map< std::string, math::Histogram> co_histograms( score::ContactOrder::HistogramsFromColumns( contact_order_table));

        // iterate over histograms
        for
        (
          storage::Map< std::string, math::Histogram>::const_iterator
            histogram_itr( co_histograms.Begin()), histogram_itr_end( co_histograms.End());
          histogram_itr != histogram_itr_end;
          ++histogram_itr
        )
        {
          // skip the first two columns
          if( histogram_itr->first == "nr_aas" || histogram_itr->first == "nr_aas_sses")
          {
            continue;
          }

          stream << histogram_itr->first + " histogram" << '\n';
          stream << histogram_itr->second << '\n';
        }
      }

      return stream.str();
    } // namespace scorestat

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ContactOrder::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes contact order statistics."
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "Histogram"
      );

      return parameters;
    }

  } // namespace scorestat
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
#include "scorestat/bcl_scorestat.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace scorestat
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
// include header for this class
#include "assemble/bcl_assemble_fold_template.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_topology_combined.h"
#include "assemble/bcl_assemble_collector_topology_interface.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "scorestat/bcl_scorestat_fold_template.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FoldTemplate::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new FoldTemplate())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &FoldTemplate::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
          "List",
          GetStaticClassName< FoldTemplate::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &FoldTemplate::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_file_extensions[] =
      {
          "fold_template.list",
          GetStaticClassName< FoldTemplate::OutputOption>()
      };
      return s_output_file_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FoldTemplate::FoldTemplate() :
        m_OutputOption( e_List),
        m_ChainIds( "")
    {
      // nothing to do
    }

    //! @brief virtual copy constructor
    FoldTemplate *FoldTemplate::Clone() const
    {
      return new FoldTemplate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &FoldTemplate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &FoldTemplate::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &FoldTemplate::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &FoldTemplate::GetAlias() const
    {
      static std::string s_name( "FoldTemplate");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string FoldTemplate::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the proein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "The protein ensemble is empty.");

      // create variables for storing output information
      storage::Vector< assemble::FoldTemplate> fold_templates;
      storage::Map
      <
        storage::Pair< size_t, size_t>,
        storage::Vector< storage::Pair< std::string, double> >
      > sorted_pdbs;

      // iterate over protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        const assemble::ProteinModel &protein_model( **protein_model_itr);

        // get pdb id
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");
        std::string pdb_id( io::File::RemovePath( model_name));
        std::transform( pdb_id.begin(), pdb_id.end(), pdb_id.begin(), toupper);

        // iterate over all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator
            chain_itr( protein_model.GetChains().Begin()), chain_itr_end( protein_model.GetChains().End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip chains that are not desired
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // create a protein model from this chain only
          assemble::ProteinModel chain_model( util::ShPtr< assemble::Chain>( ( *chain_itr)->HardCopy()));

          // create the fold template
          const util::ShPtr< assemble::CollectorTopologyInterface> sp_collector
          (
            new assemble::CollectorTopologyCombined( false)
          );
          chain_model.Transform( math::Inverse( chain_model.GetOrientation()));
          const assemble::FoldTemplate fold_template( chain_model, sp_collector, pdb_id);

          // push back the fold template information if it has at least one geometry, all defined geometries, and an appropriate Rg
          if
          (
            !fold_template.GetGeometries().IsEmpty() &&
            fold_template.HasDefinedGeometries() &&
            fold_template.GetRadiusOfGyration() < 3.0 * double( fold_template.GetGeometries().GetSize())
          )
          {
            fold_templates.PushBack( fold_template);
            sorted_pdbs
            [
               storage::Pair< size_t, size_t>
              (
                fold_template.GetHelicalGeometries().GetSize(), fold_template.GetStrandGeometries().GetSize()
              )
            ].PushBack
            (
              storage::Pair< std::string, double>( fold_template.GetPDBID(), fold_template.GetRadiusOfGyration())
            );
          }
        } // end of iterating over all chains
      } // end of iterating over the current protein model

      // write list of pdbs with fold template information
      io::OFStream write;
      io::File::MustOpenOFStream( write, "fold_template_pdbs.list");
      write << sorted_pdbs;
      io::File::CloseClearFStream( write);

      // write statistics
      std::stringstream stream;
      if( m_OutputOption == e_List)
      {
        stream << fold_templates.GetSize() << '\n';
        // iterate over fold templates
        for
        (
          storage::Vector< assemble::FoldTemplate>::const_iterator
            fold_template_itr( fold_templates.Begin()), fold_template_itr_end( fold_templates.End());
          fold_template_itr != fold_template_itr_end;
          ++fold_template_itr
        )
        {
          fold_template_itr->WriteCompact( stream);
        }
      }
      return stream.str();
    } // end of operator ()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FoldTemplate::GetSerializer() const
    {

      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes fold template statistics."
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "List"
      );
      return parameters;
    }
  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_loop_angle.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_aa.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram_2d.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "score/bcl_score_loop.h"
#include "score/bcl_score_loop_angle.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LoopAngle::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new LoopAngle())
    );

    //! @brief OutputOption as string
    //! @param OutputOption the OutputOption
    //! @return the string for the OutputOption
    const std::string &LoopAngle::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_descriptors[] =
      {
        "NormalizedByTotalCount",
        "NormalizedByColumn",
        "Log",
        "Raw",
        "Table",
        GetStaticClassName< OutputOption>()
      };

      return s_descriptors[ OUTPUT_OPTION];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &LoopAngle::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_output_file_extensions[] =
      {
        "loop_angle_count_norm.histograms",
        "loop_angle_column_norm.histograms",
        "loop_angle_log.histograms",
        "loop_angle_raw.histograms",
        "loop_angle.tbl",
        GetStaticClassName< OutputOption>()
      };

      return s_output_file_extensions[ OUTPUT_OPTION];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LoopAngle::LoopAngle() :
      m_OutputOption( e_NormalizedByCount),
      m_Chains( "A"),
      m_VisualizationFlag( false)
    {
    }

    //! @brief virtual copy constructor
    LoopAngle *LoopAngle::Clone() const
    {
      return new LoopAngle( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LoopAngle::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &LoopAngle::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LoopAngle::GetAlias() const
    {
      static const std::string s_Name( "LoopAngle");
      return s_Name;
    }

  //////////////
  // operator //
  //////////////

    //! @brief add statistics about the given protein and chains to the internal table and histograms
    //! @param ENSEMBLE the protein ensemble the statistics of which to be computed
    std::string LoopAngle::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // Initialize variables
      size_t max_sequence_distance( score::LoopAngle::GetDefaultMaxmimumSequenceDistance());

      // min value for cosine
      const double min_cos( -1.0);

      // minimum loop size
      const double min_loop_len( 0.0);

      // loop size increment
      const double loop_size_increment( 1.0);

      // Number of bins as integer
      const int x_num_bins( s_DefaultNumberBinsX);
      const int y_num_bins( s_DefaultNumberBinsY);

      // compute bin size for y axis
      const double bin_size_y( 2.0 / double( y_num_bins));

      // loops with sequence distance <= m_MaxsSequenceDistance
      math::Histogram cos_angle_short_loops_histogram( min_cos, bin_size_y, s_DefaultNumberBinsY);

      // loops with sequence distance > m_MaxsSequenceDistance
      math::Histogram cos_angle_long_loops_histogram( min_cos, bin_size_y, s_DefaultNumberBinsY);

      math::Histogram2D euclidean_distance_cos_angle_histogram
      (
        storage::VectorND< 2, double>( min_loop_len, min_cos),
        storage::VectorND< 2, double>( loop_size_increment, bin_size_y),
        storage::VectorND< 2, size_t>( x_num_bins, y_num_bins)
      );

      math::Histogram2D sequence_distance_cos_angle_histogram
      (
        storage::VectorND< 2, double>( min_loop_len, min_cos),
        storage::VectorND< 2, double>( loop_size_increment, bin_size_y),
        storage::VectorND< 2, size_t>( x_num_bins, y_num_bins)
      );

      math::Histogram2D euclidean_over_sequence_distance_cos_angle_histogram
      (
        storage::VectorND< 2, double>( min_loop_len, min_cos),
        storage::VectorND< 2, double>( loop_size_increment, bin_size_y),
        storage::VectorND< 2, size_t>( x_num_bins, y_num_bins)
      );

      storage::Table< double> loop_angle_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            "sequence_dist", "euclidean_dist", "de/log(ds)", "consecutive_sses", "cos_angle"
          )
        )
      );

      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // iterate through the ensemble
      for
      (
        assemble::ProteinEnsemble::const_iterator protein_itr( ENSEMBLE.Begin()), protein_itr_end( ENSEMBLE.End());
        protein_itr != protein_itr_end;
        ++protein_itr
      )
      {
        util::ShPtr< assemble::ProteinModel> sp_protein_model( *protein_itr);
        const assemble::ProteinModel &protein_model( *sp_protein_model);

        // get pdb filename
        const util::ShPtr< util::Wrapper< std::string> > &model_filename_ptr
        (
          protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
        );

        const std::string model_filename( model_filename_ptr.IsDefined() ? model_filename_ptr->GetData() : "");

        // create pdb factory for writing visualization pdbs
        pdb::Factory factory;

        //iterate over all chains
        const util::ShPtrVector< assemble::Chain> &chains( protein_model.GetChains());
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( chains.Begin()), chain_itr_end( chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // To simplify naming in the code
          const assemble::Chain &chain( **chain_itr);

          // skip chains that are not desired
          if( m_Chains.find( chain.GetChainID()) == std::string::npos)
          {
            BCL_MessageDbg( "Skip chain " + util::Format()( chain.GetChainID()) + std::string( ", not in chains to use: ") + m_Chains);
            continue;
          }

          // calculate radius of gyration
          const double radius_of_gyration( coord::RadiusOfGyration( chain.GetAtomCoordinates()));
          BCL_MessageDbg( "radius_of_gyration=" + util::Format()( radius_of_gyration));
          // calculate center of mass
          linal::Vector3D center_of_mass( coord::CenterOfMass( chain.GetAtomCoordinates(), true));
          BCL_MessageDbg( "center_of_mass=" + center_of_mass.ToString());

          // for each pair of (consecutive) sses: calculate euclidean distance and angles to center of gravity
          storage::Set< biol::SSType> ss_types;
          ss_types.InsertElement( biol::GetSSTypes().HELIX); // only collect helices and strands, no coil
          ss_types.InsertElement( biol::GetSSTypes().STRAND);
          util::SiPtrVector< const assemble::SSE> sses( chain.GetSSEs( ss_types));
          size_t sse_a_number( 0); // count the position of sse_a to clarify the output

          // Storage List appended atom lines
          util::ShPtrList< pdb::Line> vis_file_lines;

          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr_a( sses.Begin()), sse_itr_end( sses.End());
            sse_itr_a != sse_itr_end;
            ++sse_itr_a, ++sse_a_number
          )
          {
            if( !( *sse_itr_a)->IsDefined()) // skip sses with undefined body
            {
              continue;
            }

            // start itr for second sse one after the first
            util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr_b( sse_itr_a);
            size_t sse_b_number( sse_a_number); // count the position of sse_b to clarify the output
            bool sse_a_b_consecutive( false); // currently sse_itr_b points to the same sse as sse_itr_a
            if( sse_itr_b != sse_itr_end) // move itr to next sse if not at end
            {
              ++sse_itr_b;
              ++sse_b_number;
              sse_a_b_consecutive = true; // now it points to the next one after sse_itr_a
            }
            // iterate over all pairs of sses
            for( ; sse_itr_b != sse_itr_end; ++sse_itr_b, ++sse_b_number, sse_a_b_consecutive = false)
            {
              if( !( *sse_itr_b)->IsDefined()) // skip sses with undefined body
              {
                continue;
              }

              BCL_MessageDbg
              (
                "Calculate loop angle statistics for sse pair: sse" + util::Format()( sse_a_number) + "="
                + util::Format()( ( **sse_itr_a).GetIdentification())
                + "  sse" + util::Format()( sse_b_number) + "=" + util::Format()( ( **sse_itr_b).GetIdentification())
                + "  consecutive=" + util::Format()( sse_a_b_consecutive)
              );

              // calculate the sequence distance between the two SSEs
              const util::ShPtr< biol::AABase> last_aa_of_prev_sse( ( **sse_itr_a).GetLastAA());
              const util::ShPtr< biol::AABase> first_aa_of_next_sse( ( **sse_itr_b).GetFirstAA());
              BCL_MessageDbg( "seq_id_last_aa_of_prev_sse=" + util::Format()( last_aa_of_prev_sse->GetSeqID()));
              BCL_MessageDbg( "seq_id_first_aa_of_next_sse=" + util::Format()( first_aa_of_next_sse->GetSeqID()));

              // calculate the sequence distance
              const size_t sequence_distance( biol::SequenceSeparation( *last_aa_of_prev_sse, *first_aa_of_next_sse));
              BCL_MessageDbg( "sequence_distance=" + util::Format()( sequence_distance));

              // get positions of begin and end of z-axis of the two SSEs
              const linal::Vector3D prev_end_of_z( ( **sse_itr_a).EndOfZ());
              const linal::Vector3D next_begin_of_z( ( **sse_itr_b).BeginOfZ());
              BCL_MessageDbg( "prev_end_of_z=" + prev_end_of_z.ToString());
              BCL_MessageDbg( "next_begin_of_z=" + next_begin_of_z.ToString());

              // calculate the euclidean distance
              const double euclidean_distance( linal::Distance( next_begin_of_z, prev_end_of_z));
              BCL_MessageDbg( "euclidean_distance=" + util::Format()( euclidean_distance));

              // cosine of projection angle between center_of_mass->prev_end_of_z and center_of_mass->next_begin_of_z
              const double cosine_of_proj_angle( linal::ProjAngleCosinus( center_of_mass, prev_end_of_z, next_begin_of_z));
              BCL_MessageDbg( "cosine_of_proj_angle=" + util::Format()( cosine_of_proj_angle));

              if( m_OutputOption != e_Table)
              {
                // add data to histograms
                if( sequence_distance <= max_sequence_distance)
                {
                  cos_angle_short_loops_histogram.PushBack( cosine_of_proj_angle);
                }
                else
                {
                  cos_angle_long_loops_histogram.PushBack( cosine_of_proj_angle);
                }

                euclidean_distance_cos_angle_histogram.PushBack
                (
                  storage::VectorND< 2, double>( euclidean_distance, cosine_of_proj_angle)
                );

                sequence_distance_cos_angle_histogram.PushBack
                (
                  storage::VectorND< 2, double>( sequence_distance, cosine_of_proj_angle)
                );

                euclidean_over_sequence_distance_cos_angle_histogram.PushBack
                (
                  storage::VectorND< 2, double>( euclidean_distance / std::log( sequence_distance), cosine_of_proj_angle)
                );
              }
              else
              {
                // add data to table
                const std::string sse_str( "sse" + util::Format()( sse_a_number) + "_sse" + util::Format()( sse_b_number));
                loop_angle_table.InsertRow
                (
                  model_filename + "_" + sse_str,
                  storage::Vector< double>::Create
                  (
                    sequence_distance,
                    euclidean_distance,
                    score::Loop::NormalizeDistance( storage::Pair< size_t, double>( sequence_distance, euclidean_distance)),
                    sse_a_b_consecutive,
                    cosine_of_proj_angle
                  ),
                  true
                );
              }

              // only perform if path was set from commandline
              if( m_VisualizationFlag)
              {
                // add atom lines for center_of_mass, prev_end_of_z, next_begin_of_z
                biol::Atom atom_center_of_mass( center_of_mass, biol::GetAtomTypes().CA);
                biol::Atom atom_prev_end_of_z( prev_end_of_z, biol::GetAtomTypes().CA);
                biol::Atom atom_next_begin_of_z( next_begin_of_z, biol::GetAtomTypes().CA);
                biol::AA amino_acid;

                vis_file_lines.Append
                (
                  pdb::Factory::WriteAtomToLine( atom_center_of_mass, amino_acid, 'Z', 1 + vis_file_lines.GetSize())
                );

                vis_file_lines.Append
                (
                  pdb::Factory::WriteAtomToLine( atom_prev_end_of_z, amino_acid, 'Z', 1 + vis_file_lines.GetSize())
                );

                vis_file_lines.Append
                (
                  pdb::Factory::WriteAtomToLine( atom_next_begin_of_z, amino_acid, 'Z', 1 + vis_file_lines.GetSize())
                );

              } // end if

            } // for sse_itr_b
          } // for sse_itr_a

          // Write the file for each protein in the ensemble
          if( m_VisualizationFlag)
          {
            std::string vis_filename( model_filename + GetOutFilePostfix() + ".pdb");

            // create handler and add lines
            pdb::Handler pdb_handler;
            pdb_handler.AppendLines( vis_file_lines);

            // write visualization pdb
            io::OFStream pdb_write_stream;
            io::File::MustOpenOFStream( pdb_write_stream, vis_filename);
            pdb_handler.WriteLines( pdb_write_stream);
            io::File::CloseClearFStream( pdb_write_stream);
          }

        } // for chain
      }// End Ensemble iteration

      std::stringstream ostream;
      if( m_OutputOption == e_Table)
      {
        loop_angle_table.WriteFormatted( ostream);
      }
      else
      {
        // normalization for histogram outputs
        switch( m_OutputOption)
        {
          case e_NormalizedByCount:
            euclidean_distance_cos_angle_histogram.Normalize();
            sequence_distance_cos_angle_histogram.Normalize();
            euclidean_over_sequence_distance_cos_angle_histogram.Normalize();
            break;
          case e_NormalizedByColumn:
            euclidean_distance_cos_angle_histogram.NormalizeY();
            sequence_distance_cos_angle_histogram.NormalizeY();
            euclidean_over_sequence_distance_cos_angle_histogram.NormalizeY();
            break;
          case e_Log:
            euclidean_distance_cos_angle_histogram.Log();
            sequence_distance_cos_angle_histogram.Log();
            euclidean_over_sequence_distance_cos_angle_histogram.Log();
            break;
          default:
            break;
        } // end normalization
        ostream << cos_angle_short_loops_histogram << '\n'
                << cos_angle_long_loops_histogram << '\n'
                << euclidean_distance_cos_angle_histogram << '\n'
                << sequence_distance_cos_angle_histogram << '\n'
                << euclidean_over_sequence_distance_cos_angle_histogram;
      }
      return ostream.str();
    } // end of operator()

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LoopAngle::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes loop distance statistics."
      );

      parameters.AddInitializer
      (
        "output option",
        "the form of output to create during analysis, Normalized by count, Normalized by Column, Log, or Raw",
        io::Serialization::GetAgent( &m_OutputOption),
        "NormalizedByTotalCount"
      );

      parameters.AddInitializer
      (
        "chains",
        "a string of chains to use for the analysis",
        io::Serialization::GetAgent( &m_Chains),
        "A"
      );

      parameters.AddInitializer
      (
        "visualize",
        "set path to visualize the model in pymol",
        io::Serialization::GetAgent( &m_VisualizationFlag),
        "False"
      );

      return parameters;
    } // end of GetSerializer function
  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_loop_closure.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LoopClosure::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new LoopClosure())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &LoopClosure::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Histogram",
        GetStaticClassName< LoopClosure::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &LoopClosure::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "loop_closure.histogram",
        GetStaticClassName< LoopClosure::OutputOption>()
      };
      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LoopClosure::LoopClosure() :
      m_OutputOption( e_Histogram),
      m_DistanceBinSize( 0.1),
      m_MaxDistance( 50.001), //
      m_NumResidues( 30),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    LoopClosure *LoopClosure::Clone() const
    {
      return new LoopClosure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &LoopClosure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &LoopClosure::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns the bin size for the histogram
    //! @return the bin size for the histogram
    const double &LoopClosure::GetDistanceBinSize() const
    {
      return m_DistanceBinSize;
    }

    //! @brief returns the maximum distance
    //! @returns the maximum distance
    const double &LoopClosure::GetMaxDistance() const
    {
      return m_MaxDistance;
    }

    //! @brief returns the number of residues
    //! @returns the number of residues
    const size_t &LoopClosure::GetNumResidues() const
    {
      return m_NumResidues;
    }

    //! @brief returns chain id
    //! @return chain id
    const std::string &LoopClosure::GetChainId() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LoopClosure::GetAlias() const
    {
      static std::string s_name( "LoopClosure");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string LoopClosure::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // statistics for loop closure
      storage::Vector< math::Histogram> loop_closure_histograms
      (
        m_NumResidues,
        math::Histogram( 0.0, m_DistanceBinSize, size_t( m_MaxDistance / m_DistanceBinSize))
      );

      // iterator over the protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get the current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get protein pdb filename
        const util::ShPtr< util::Wrapper< std::string> >
          &sp_model_filename( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_filename( sp_model_filename.IsDefined() ? io::File::RemovePath( sp_model_filename->GetData()) : "");

        // get all chains in current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // iterate over all chains in current protein model
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator
            chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // get all sses in current chain
          const storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> &all_sses( ( *chain_itr)->GetData());

          // iterate over all sses in current chain
          for
          (
            storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
              sse_itr_a( all_sses.Begin()), sse_itr_end( all_sses.End());
            sse_itr_a != sse_itr_end;
            ++sse_itr_a
          )
          {
            // get current sse type
            const biol::SSType &current_sse_type( ( *sse_itr_a)->GetType());

            // skip protein models without sse entries in the pdb file
            if( ( *chain_itr)->GetNumberSSEs() == 1 && current_sse_type == biol::GetSSTypes().COIL)
            {
              BCL_MessageStd( util::Format()( model_filename) + " probably does not have sse entries in the pdb file, skipping");
              continue;
            }

            // skip undefined sses
            if( !( *sse_itr_a)->IsDefined())
            {
              continue;
            }

            // iterate over the second sse
            storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator sse_itr_b( sse_itr_a);
            ++sse_itr_b;
            for
            (
              ;
              sse_itr_b != sse_itr_end;
              ++sse_itr_b
            )
            {
              // skip undefined sses
              if( !( *sse_itr_b)->IsDefined())
              {
                continue;
              }

              // get the sequence distance
              const size_t sequence_distance( biol::CalculateSequenceDistance( **sse_itr_a, **sse_itr_b));

              // if the distance is larger than needed
              if( sequence_distance > loop_closure_histograms.GetSize() - 1)
              {
                continue;
              }

              // get euclidean distance
              const double euclidean_distance
              (
                biol::Distance
                (
                  ( *sse_itr_a)->GetLastAA()->GetAtom( biol::GetAtomTypes().C),
                  ( *sse_itr_b)->GetFirstAA()->GetAtom( biol::GetAtomTypes().N)
                )
              );

              // skip undefined distance
              if( euclidean_distance == util::GetUndefinedDouble())
              {
                continue;
              }

              // push the distance into histogram
              if( m_OutputOption == e_Histogram)
              {
                loop_closure_histograms( sequence_distance).PushBack( euclidean_distance);
              }

            } // end of iteration over sses
          } // end of iteration over chains
        } // end of iteration over current protein model
      } // end of iteration over protein ensemble

      // write statistics
      std::stringstream stream;

      if( m_OutputOption == e_Histogram)
      {
        stream << loop_closure_histograms;
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LoopClosure::GetSerializer() const
    {
       io::Serializer parameters;
       parameters.SetClassDescription
       (
         "Computes loop closure statistics."
       );

       parameters.AddInitializer
       (
         "bin_size",
         "the bin size for the histogram",
         io::Serialization::GetAgent( &m_DistanceBinSize),
         "0.1"
       );

       parameters.AddInitializer
       (
         "maximum_distance",
         "maximum distance",
         io::Serialization::GetAgent( &m_MaxDistance),
         "50.001"
       );

       parameters.AddInitializer
       (
         "number_residues",
         "number of residues",
         io::Serialization::GetAgent( &m_NumResidues),
         "30"
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
         "output",
         "what type of outputs to provide",
         io::Serialization::GetAgent( &m_OutputOption),
         "Histogram"
       );

       return parameters;
    }
  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_loop_distance.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score_loop.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LoopDistance::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new LoopDistance())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &LoopDistance::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_Names[] =
      {
        "Table",
        "Histogram",
        "LogHistogram",
        GetStaticClassName< LoopDistance::OutputOption>()
      };
      return s_Names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &LoopDistance::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_output_file_extensions[] =
      {
        "loop_distance.tbl",
        "loop_distance_raw.histograms",
        "loop_distance_log.histograms",
        GetStaticClassName< OutputOption>()
      };

      return s_output_file_extensions[ OUTPUT_OPTION];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LoopDistance::LoopDistance() :
      m_OutputOption( e_Table),
      m_OutFilePostFix( ".loopdistance"),
      m_BinSize( 1.0),
      m_ChainIds( "")
    {
    }

    //! @brief virtual copy constructor
    LoopDistance *LoopDistance::Clone() const
    {
      return new LoopDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LoopDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &LoopDistance::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns the binsize for the histogram
    //! @return the binsize for the histogram
    const double &LoopDistance::GetBinSize() const
    {
      return m_BinSize;
    }

    //! @brief returns chain id
    //! @return chain id
    const std::string &LoopDistance::GetChainId() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LoopDistance::GetAlias() const
    {
      static const std::string s_Name( "LoopDistance");
      return s_Name;
    }

  //////////////
  // operator //
  //////////////

    //! @brief add statistics about the given protein and chains to the internal table and histograms
    //! @param ENSEMBLE the protein ensemble the statistics of which to be computed
    std::string LoopDistance::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // three output formats
      storage::Vector< math::Histogram> loop_dist_histogram
      (
        500, math::Histogram( double( 0), m_BinSize, size_t( 100 / m_BinSize))
      );

      storage::Table< double> loop_dist_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            "sequence_dist", "euclidean_dist", "de/log(ds)", "consecutive_sses"
          )
        )
      );

      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // iterate over the protein ensemble
      for
      (
         util::ShPtrVector< assemble::ProteinModel>::const_iterator
           protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
         protein_model_itr != protein_model_itr_end;
         ++protein_model_itr
      )
      {
        // get current protein model
        util::ShPtr< assemble::ProteinModel> sp_protein_model( *protein_model_itr);
        const assemble::ProteinModel &protein_model( *sp_protein_model);

        // get current pdb name before all the loops start
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string model_name( model_name_ptr.IsDefined() ? io::File::RemovePath( model_name_ptr->GetData()) : "");

        // iterate over all chains of the protein model
        for
        (
           util::ShPtrVector< assemble::Chain>::const_iterator
             chain_itr( protein_model.GetChains().Begin()), chain_itr_end( protein_model.GetChains().End());
           chain_itr != chain_itr_end;
           ++chain_itr
        )
        {
          // get current chain
          util::ShPtr< assemble::Chain> current_chain( *chain_itr);

          BCL_MessageDbg
          (
            "Calculating loop distance statistics for chain, chain.sse.size="
            + util::Format()( current_chain->GetData().GetSize())
          );

          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( current_chain->GetChainID()) == std::string::npos)
          {
            BCL_MessageDbg
            (
              std::string( "Skip chain: ") + current_chain->GetChainID() + std::string( ", chain ids to use: ") + m_ChainIds
            );

            continue;
          }

          // iterate over all sses in the current chain
          size_t sse_a_number( 0);
          for
          (
             storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
               sse_a_itr( current_chain->GetData().Begin()), sse_itr_end( current_chain->GetData().End());
              sse_a_itr != sse_itr_end;
              ++sse_a_itr
          )
          {

            // get the first sse of the current chain
            storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
              sse_b_itr( sse_a_itr);

            // skip undefined sse
            if( !( *sse_a_itr)->IsDefined())
            {
              continue;
            }

            size_t sse_b_number( sse_a_number);
            bool is_sse_a_b_consecutive( false);
            if( sse_b_itr != sse_itr_end)
            {
              ++sse_b_itr;
              ++sse_b_number;
              is_sse_a_b_consecutive = true;
            }

            // iterate over the rest sses
            for( ; sse_b_itr != sse_itr_end; ++sse_b_itr)
            {
              // skip undefined sse
              if( !( *sse_b_itr)->IsDefined())
              {
                BCL_MessageDbg( "Skip: sse_b is undefined");
                continue;
              }

              BCL_MessageDbg
              (
                "Calculate loop distance statistics for sse pair " +
                util::Format()( ( *sse_a_itr)->GetIdentification()) +
                " " + util::Format()( ( *sse_b_itr)->GetIdentification())
              );

              // calculate sequence distance and euclidean distance
              const storage::Pair< size_t, double> seq_euc_distance
              (
                score::Loop::SequenceAndEuclideanDistance( **sse_a_itr, **sse_b_itr)
              );

              // store sequence distance and euclidean distance in histogram
              if( seq_euc_distance.First() >= loop_dist_histogram.GetSize())
              {
                if( m_OutputOption == e_Histogram)
                {
                  loop_dist_histogram.LastElement().PushBack( seq_euc_distance.Second());
                }
                else if( m_OutputOption == e_LogHistogram)
                {
                  loop_dist_histogram.LastElement().PushBack( score::Loop::NormalizeDistance( seq_euc_distance));
                }
              }
              else
              {
                if( m_OutputOption == e_Histogram)
                {
                  loop_dist_histogram( seq_euc_distance.First()).PushBack( seq_euc_distance.Second());
                }
                else if( m_OutputOption == e_LogHistogram)
                {
                  loop_dist_histogram( seq_euc_distance.First()).PushBack( score::Loop::NormalizeDistance( seq_euc_distance));
                }
              }

              // add loop distance statistics to table
              if
              (
                  ( **sse_a_itr).GetType() != biol::GetSSTypes().COIL
                  && ( **sse_b_itr).GetType() != biol::GetSSTypes().COIL
              )
              {
                BCL_MessageDbg
                (
                  "Insert into table: seq_dist=" + util::Format()( seq_euc_distance.First())
                  + "  euc_dist=" + util::Format()( seq_euc_distance.Second())
                  + "  norm_dist=" + util::Format()( score::Loop::NormalizeDistance( seq_euc_distance))
                );

                if( m_OutputOption == e_Table)
                {
                  const std::string sse_pair( "sse" + util::Format()( sse_a_number) + "_sse" + util::Format()( sse_b_number));
                  loop_dist_table.InsertRow
                  (
                    model_name + "_" + sse_pair,
                    storage::Vector< double>::Create
                    (
                      seq_euc_distance.First(),
                      seq_euc_distance.Second(),
                      score::Loop::NormalizeDistance( seq_euc_distance),
                      is_sse_a_b_consecutive
                    ),
                    true
                  );
                }

                BCL_MessageDbg( "table size=" + util::Format()( loop_dist_table.GetSize()));

                is_sse_a_b_consecutive = false;
                ++sse_b_number;
              }
            }

            if( ( **sse_a_itr).GetType() != biol::GetSSTypes().COIL)
            {
              ++sse_a_number;
            }
          }
        } // end of the current chain
      } // end of current model

      // write statistics
      std::ostringstream stream;
      if( m_OutputOption == e_Histogram || m_OutputOption == e_LogHistogram)
      {
        stream << loop_dist_histogram;

      }
      else if( m_OutputOption == e_Table)
      {
        loop_dist_table.WriteFormatted( stream);
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LoopDistance::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes loop distance statistics."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".loopdistance"
      );

      parameters.AddInitializer
      (
        "bin_size",
        "the bin size for the histogram",
        io::Serialization::GetAgent( &m_BinSize)
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "Table"
      );

      return parameters;
    }

  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_neighbor_count.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> NeighborCount::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new NeighborCount( false))
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> NeighborCount::s_InstanceEnv
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new NeighborCount( true))
    );

    //! @brief ChainOption as string
    //! @param CHAIN_OPTION the ChainOption
    //! @return the string for the ChainOption
    const std::string &NeighborCount::GetChainOptionName( const ChainOption &CHAIN_OPTION)
    {
      static std::string s_names[] =
      {
        "One",
        "All",
        GetStaticClassName< NeighborCount::ChainOption>()
      };
      return s_names[ size_t( CHAIN_OPTION)];
    }

    //! @brief Output filename as string
    //! @param ChainOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &NeighborCount::GetOutputFileName( const ChainOption &CHAIN_OPTION, const bool SPLIT_ENVIRONMENT)
    {
      static std::string s_names_default[] =
      {
        "nc_one_chain.wo_env.histograms",
        "nc_all_chain.wo_env.histograms",
        GetStaticClassName< NeighborCount::ChainOption>()
      };

      static std::string s_names_env[] =
      {
        "nc_one_chain.with_env.histograms",
        "nc_all_chain.with_env.histograms",
        GetStaticClassName< NeighborCount::ChainOption>()
      };

      return SPLIT_ENVIRONMENT ? s_names_env[ size_t( CHAIN_OPTION)] : s_names_default[ size_t( CHAIN_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    NeighborCount::NeighborCount()
    {
      // nothing else to do
    }

    //! @brief explicit constructor
    NeighborCount::NeighborCount( const bool SPLIT_ENVIRONMENT) :
        m_ChainOption( e_OneChain),
        m_ChainIds( ""),
        m_SequenceExclusion( 2),
        m_NCLowerBound( 4.0),
        m_NCUpperBound( 11.4),
        m_SplitEnvironment( SPLIT_ENVIRONMENT)
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    NeighborCount *NeighborCount::Clone() const
    {
      return new NeighborCount( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &NeighborCount::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &NeighborCount::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_ChainOption, m_SplitEnvironment);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &NeighborCount::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the sequence exclusion considered when creating neighbor list
    //! @return the sequence exclusion considered when creating neighbor list
    const size_t &NeighborCount::GetSequenceExclusion() const
    {
      return m_SequenceExclusion;
    }

    //! @brief returns the lower bound for neighbor count
    //! @return the lower bound for neighbor count
    const double &NeighborCount::GetNCLowerBound() const
    {
      return m_NCLowerBound;
    }

    //! @brief returns the upper bound for neighbor count
    //! @return the upper bound for neighbor count
    const double &NeighborCount::GetNCUpperBound() const
    {
      return m_NCUpperBound;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &NeighborCount::GetAlias() const
    {
      static std::string s_name_one( "NeighborCount"), s_name_two( "NeighborCountByEnvironment");
      return m_SplitEnvironment ? s_name_two : s_name_one;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string NeighborCount::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // statistics for neighbor count one chain, categorized by environment types
      storage::Vector< storage::Vector< math::Histogram> > desired_histograms
      (
        biol::GetEnvironmentTypes().GetEnumCount(),
        storage::Vector< math::Histogram>( biol::GetAATypes().GetEnumCount(), math::Histogram( 0, 1.0, 50))
      );

      // initialize a neighbor count object
      assemble::AANeighborCount neighbor_count;
      neighbor_count.SetMinimalSequenceSeparation( m_SequenceExclusion);
      neighbor_count.SetThresholdRange( math::Range< double>( m_NCLowerBound, m_NCUpperBound));

      // iterate over the protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get the membrane
        const util::ShPtr< biol::Membrane> sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        // get all chains of current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // get all amino acids and cb coordinates
        util::SiPtrVector< const biol::AABase> all_chain_amino_acids;

        // iterate over all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
            chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // get all amino acids
          all_chain_amino_acids.Append( ( *chain_itr)->GetAminoAcids());
        }

        // get neighbor lists from all chains for computing neighbor vector statistics
        const assemble::AANeighborListContainer all_chain_neighbor_list_nc
        (
          all_chain_amino_acids,
          neighbor_count.GetDistanceCutoff(),
          neighbor_count.GetMinimalSequenceSeparation(),
          true
        );

        // iterate through all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
            chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // ship undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // get amino acids in current chain
          const util::SiPtrVector< const biol::AABase> amino_acids_current_chain( ( *chain_itr)->GetAminoAcids());

          // get neighbor lists from current chain for computing neighbor vector statistics
          const assemble::AANeighborListContainer current_chain_neighbor_list_nc
          (
            amino_acids_current_chain,
            neighbor_count.GetDistanceCutoff(),
            m_SequenceExclusion,
            false
          );

          const assemble::AANeighborListContainer &desired_neighbor_list
          (
            m_ChainOption == e_OneChain
            ? current_chain_neighbor_list_nc
            : all_chain_neighbor_list_nc
          );

          // iterate through current chain
          for
          (
            util::SiPtrVector< const biol::AABase>::const_iterator
              aa_itr( amino_acids_current_chain.Begin()), aa_itr_end( amino_acids_current_chain.End());
            aa_itr != aa_itr_end;
            ++aa_itr
          )
          {
            // proceed only if there are coordinates given and the aa is a natural aa
            if
            (
              !( *aa_itr)->GetFirstSidechainAtom().GetCoordinates().IsDefined() ||
              !( *aa_itr)->GetType()->IsNaturalAminoAcid()
            )
            {
              continue;
            }

            // get current environment type
            const biol::EnvironmentType current_environment_type
            (
              sp_membrane.IsDefined() && m_SplitEnvironment ?
                  sp_membrane->DetermineEnvironmentType( ( *aa_itr)->GetFirstSidechainAtom().GetCoordinates()) :
                  biol::GetEnvironmentTypes().e_Solution
            );

            // calculate desired neighbor count
            const double &desired_neighbor_count
            (
              neighbor_count( desired_neighbor_list.Find( ( **aa_itr))->second)
            );

            desired_histograms( current_environment_type)( ( *aa_itr)->GetType()).PushBack( desired_neighbor_count);

          } // end of iterating through all aas in current chain
        } // end of iterating through all chains
      } // end of iterating through all protein models

      // write output
      std::stringstream stream;

      stream << neighbor_count.GetThresholdRange() << '\n';
      stream << neighbor_count.GetMinimalSequenceSeparation() << '\n';

      // iterate over environment types
      for
      (
        biol::EnvironmentTypes::const_iterator
          ent_itr( biol::GetEnvironmentTypes().Begin()), ent_itr_end( biol::GetEnvironmentTypes().End());
        ent_itr != ent_itr_end;
        ++ent_itr
      )
      {
        // only interested in soluble environment type
        if( !m_SplitEnvironment)
        {
          if( ( *ent_itr) != biol::GetEnvironmentTypes().e_Solution)
          {
            continue;
          }
        }
        else
        {
          stream << *ent_itr << '\n';
        }

        // for current environment type iterate over all amino acid types
        for
        (
          biol::AATypes::const_iterator
            aa_itr( biol::GetAATypes().Begin()),
            aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          stream << ( *aa_itr)->GetOneLetterCode() << '\n';
          stream << desired_histograms( *ent_itr)( *aa_itr);
        }
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NeighborCount::GetSerializer() const
    {

      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes neighbor count statistics."
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "sequence_exclusion",
        "minimal sequence separation to be considered in the calculation",
        io::Serialization::GetAgent( &m_SequenceExclusion),
        "2"
      );

      parameters.AddInitializer
      (
        "lower_bound",
        "lower bound for distance considered in neighbor count calculation",
        io::Serialization::GetAgent( &m_NCLowerBound),
        "4.0"
      );

      parameters.AddInitializer
      (
        "upper_bound",
        "upper bound for distance considered in neighbor count calculation",
        io::Serialization::GetAgent( &m_NCUpperBound),
        "11.4"
      );

      parameters.AddInitializer
      (
        "chain_option",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_ChainOption),
        GetChainOptionName( e_OneChain)
      );
      return parameters;
    }
  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_neighbor_vector.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> NeighborVector::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new NeighborVector( false))
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> NeighborVector::s_InstanceEnv
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new NeighborVector( true))
    );

    //! @brief ChainOption as string
    //! @param CHAIN_OPTION the ChainOption
    //! @return the string for the ChainOption
    const std::string &NeighborVector::GetChainOptionName( const ChainOption &CHAIN_OPTION)
    {
      static std::string s_names[] =
      {
        "One",
        "All",
        GetStaticClassName< NeighborVector::ChainOption>()
      };
      return s_names[ size_t( CHAIN_OPTION)];
    }

    //! @brief Output filename as string
    //! @param ChainOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &NeighborVector::GetOutputFileName( const ChainOption &CHAIN_OPTION, const bool SPLIT_ENVIRONMENT)
    {
      static std::string s_names_default[] =
      {
        "nv_one_chain.wo_env.histograms",
        "nv_all_chain.wo_env.histograms",
        GetStaticClassName< NeighborVector::ChainOption>()
      };

      static std::string s_names_env[] =
      {
        "nv_one_chain.with_env.histograms",
        "nv_all_chain.with_env.histograms",
        GetStaticClassName< NeighborVector::ChainOption>()
      };

      return SPLIT_ENVIRONMENT ? s_names_env[ size_t( CHAIN_OPTION)] : s_names_default[ size_t( CHAIN_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    NeighborVector::NeighborVector()
    {
      // nothing else to do
    }

    //! @brief explicit constructor
    NeighborVector::NeighborVector( const bool SPLIT_ENVIRONMENT) :
        m_ChainOption( e_OneChain),
        m_ChainIds( ""),
        m_SequenceExclusion( 2),
        m_NVLowerBound( 3.3),
        m_NVUpperBound( 11.1),
        m_SplitEnvironment( SPLIT_ENVIRONMENT)
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    NeighborVector *NeighborVector::Clone() const
    {
      return new NeighborVector( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &NeighborVector::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &NeighborVector::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_ChainOption, m_SplitEnvironment);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &NeighborVector::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the sequence exclusion considered when creating neighbor list
    //! @return the sequence exclusion considered when creating neighbor list
    const size_t &NeighborVector::GetSequenceExclusion() const
    {
      return m_SequenceExclusion;
    }

    //! @brief returns the lower bound for neighbor count
    //! @return the lower bound for neighbor count
    const double &NeighborVector::GetNVLowerBound() const
    {
      return m_NVLowerBound;
    }

    //! @brief returns the upper bound for neighbor count
    //! @return the upper bound for neighbor count
    const double &NeighborVector::GetNVUpperBound() const
    {
      return m_NVUpperBound;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &NeighborVector::GetAlias() const
    {
      static std::string s_name_one( "NeighborVector"), s_name_two( "NeighborVectorByEnvironment");
      return m_SplitEnvironment ? s_name_two : s_name_one;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string NeighborVector::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // statistics for neighbor vector one chain, categorized by environment types
      storage::Vector< storage::Vector< math::Histogram> > desired_histograms
      (
        biol::GetEnvironmentTypes().GetEnumCount(),
        storage::Vector< math::Histogram>( biol::GetAATypes().GetEnumCount(), math::Histogram( 0, 0.02, 50))
      );

      // initialize a neighbor vector object
      assemble::AANeighborVector neighbor_vector;
      neighbor_vector.SetMinimalSequenceSeparation( m_SequenceExclusion);
      neighbor_vector.SetThresholdRange( math::Range< double>( m_NVLowerBound, m_NVUpperBound));

      // iterate over the protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get the membrane
        const util::ShPtr< biol::Membrane> sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        // get all chains of current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // get all amino acids and cb coordinates
        util::SiPtrVector< const biol::AABase> all_chain_amino_acids;

        // iterate over all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
            chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // get all amino acids
          all_chain_amino_acids.Append( ( *chain_itr)->GetAminoAcids());
        }

        // get neighbor lists from all chains for computing neighbor vector statistics
        const assemble::AANeighborListContainer all_chain_neighbor_list_nc
        (
          all_chain_amino_acids,
          neighbor_vector.GetDistanceCutoff(),
          neighbor_vector.GetMinimalSequenceSeparation(),
          true
        );

        // iterate through all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
            chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // ship undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // get amino acids in current chain
          const util::SiPtrVector< const biol::AABase> amino_acids_current_chain( ( *chain_itr)->GetAminoAcids());

          // get neighbor lists from current chain for computing neighbor vector statistics
          const assemble::AANeighborListContainer current_chain_neighbor_list_nc
          (
            amino_acids_current_chain,
            neighbor_vector.GetDistanceCutoff(),
            neighbor_vector.GetMinimalSequenceSeparation(),
            false
          );

          const assemble::AANeighborListContainer &desired_neighbor_list
          (
            m_ChainOption == e_OneChain
            ? current_chain_neighbor_list_nc
            : all_chain_neighbor_list_nc
          );

          // iterate through current chain
          for
          (
            util::SiPtrVector< const biol::AABase>::const_iterator
              aa_itr( amino_acids_current_chain.Begin()), aa_itr_end( amino_acids_current_chain.End());
            aa_itr != aa_itr_end;
            ++aa_itr
          )
          {
            // proceed only if there are coordinates given and the aa is a natural aa
            if
            (
              !( *aa_itr)->GetFirstSidechainAtom().GetCoordinates().IsDefined() ||
              !( *aa_itr)->GetType()->IsNaturalAminoAcid()
            )
            {
              continue;
            }

            // get current environment type
            const biol::EnvironmentType current_environment_type
            (
              sp_membrane.IsDefined() && m_SplitEnvironment ?
                  sp_membrane->DetermineEnvironmentType( ( *aa_itr)->GetFirstSidechainAtom().GetCoordinates()) :
                  biol::GetEnvironmentTypes().e_Solution
            );

            // calculate desired neighbor vector
            const double &desired_neighbor_vector
            (
              neighbor_vector( desired_neighbor_list.Find( ( **aa_itr))->second)
            );

            desired_histograms( current_environment_type)( ( *aa_itr)->GetType()).PushBack( desired_neighbor_vector);

          } // end of iterating through all aas in current chain
        } // end of iterating through all chains
      } // end of iterating through all protein models

      // write output
      std::stringstream stream;

      stream << neighbor_vector.GetThresholdRange() << '\n';
      stream << neighbor_vector.GetMinimalSequenceSeparation() << '\n';

      // iterate over environment types
      for
      (
        biol::EnvironmentTypes::const_iterator
          ent_itr( biol::GetEnvironmentTypes().Begin()), ent_itr_end( biol::GetEnvironmentTypes().End());
        ent_itr != ent_itr_end;
        ++ent_itr
      )
      {
        // only interested in soluble environment type
        if( !m_SplitEnvironment)
        {
          if( ( *ent_itr) != biol::GetEnvironmentTypes().e_Solution)
          {
            continue;
          }
        }
        else
        {
          stream << *ent_itr << '\n';
        }

        // for current environment type iterate over all amino acid types
        for
        (
          biol::AATypes::const_iterator
            aa_itr( biol::GetAATypes().Begin()),
            aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          stream << ( *aa_itr)->GetOneLetterCode() << '\n';
          stream << desired_histograms( *ent_itr)( *aa_itr);
        }
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NeighborVector::GetSerializer() const
    {

      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes neighbor vector statistics."
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "sequence_exclusion",
        "minimal sequence separation to be considered in the calculation",
        io::Serialization::GetAgent( &m_SequenceExclusion),
        "2"
      );

      parameters.AddInitializer
      (
        "lower_bound",
        "lower bound for distance considered in neighbor count calculation",
        io::Serialization::GetAgent( &m_NVLowerBound),
        "3.3"
      );

      parameters.AddInitializer
      (
        "upper_bound",
        "upper bound for distance considered in neighbor count calculation",
        io::Serialization::GetAgent( &m_NVUpperBound),
        "11.1"
      );

      parameters.AddInitializer
      (
        "chain_option",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_ChainOption),
        GetChainOptionName( e_OneChain)
      );
      return parameters;
    }
  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_ols.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_sasa_ols.h"
#include "assemble/bcl_assemble_chain.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "math/bcl_math_histogram.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> OLS::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new OLS( false))
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> OLS::s_InstanceEnv
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new OLS( true))
    );

    //! @brief ChainOption as string
    //! @param CHAIN_OPTION the ChainOption
    //! @return the string for the ChainOption
    const std::string &OLS::GetChainOptionName( const ChainOption &CHAIN_OPTION)
    {
      static std::string s_names[] =
      {
        "One",
        "All",
        GetStaticClassName< OLS::ChainOption>()
      };
      return s_names[ size_t( CHAIN_OPTION)];
    }

    //! @brief Output filename as string
    //! @param ChainOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &OLS::GetOutputFileName( const ChainOption &CHAIN_OPTION, const bool SPLIT_ENVIRONMENT)
    {
      static std::string s_names_default[] =
      {
        "ols_one_chain.wo_env.histograms",
        "ols_all_chain.wo_env.histograms",
        GetStaticClassName< OLS::ChainOption>()
      };

      static std::string s_names_env[] =
      {
        "ols_one_chain.with_env.histograms",
        "ols_all_chain.with_env.histograms",
        GetStaticClassName< OLS::ChainOption>()
      };

      return SPLIT_ENVIRONMENT ? s_names_env[ size_t( CHAIN_OPTION)] : s_names_default[ size_t( CHAIN_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    OLS::OLS()
    {
      // nothing else to do
    }

    //! @brief explicit constructor
    OLS::OLS( const bool SPLIT_ENVIRONMENT) :
        m_ChainOption( e_OneChain),
        m_ChainIds( ""),
        m_SequenceExclusion( 2),
        m_Radius( 4.75),
        m_SplitEnvironment( SPLIT_ENVIRONMENT)
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    OLS *OLS::Clone() const
    {
      return new OLS( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &OLS::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &OLS::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_ChainOption, m_SplitEnvironment);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &OLS::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the sequence exclusion considered when creating neighbor list
    //! @return the sequence exclusion considered when creating neighbor list
    const size_t &OLS::GetSequenceExclusion() const
    {
      return m_SequenceExclusion;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &OLS::GetAlias() const
    {
      static std::string s_name_one( "OLS"), s_name_two( "OLSByEnvironment");
      return m_SplitEnvironment ? s_name_two : s_name_one;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string OLS::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // statistics for neighbor count one chain, categorized by environment types
      storage::Vector< storage::Vector< math::Histogram> > desired_histograms
      (
        biol::GetEnvironmentTypes().GetEnumCount(),
        storage::Vector< math::Histogram>( biol::GetAATypes().GetEnumCount(), math::Histogram( 0, 0.01, 100))
      );

      // initialize an ols object
      assemble::AASasaOLS ols;
      ols.SetMinimalSequenceSeparation( m_SequenceExclusion);
      ols.SetThresholdRange( math::Range< double>( 0, m_Radius));

      // iterate over the protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get the membrane
        const util::ShPtr< biol::Membrane> sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        // get all chains of current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // get all amino acids and cb coordinates
        util::SiPtrVector< const biol::AABase> all_chain_amino_acids;

        // iterate over all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
            chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // get all amino acids
          all_chain_amino_acids.Append( ( *chain_itr)->GetAminoAcids());
        }

        // get neighbor lists from all chains for computing neighbor vector statistics
        const assemble::AANeighborListContainer all_chain_neighbor_list_container
        (
          all_chain_amino_acids,
          ols.GetDistanceCutoff(),
          ols.GetMinimalSequenceSeparation(),
          true
        );

        // iterate through all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
            chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // ship undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // get amino acids in current chain
          const util::SiPtrVector< const biol::AABase> amino_acids_current_chain( ( *chain_itr)->GetAminoAcids());

          // get neighbor lists from current chain for computing neighbor vector statistics
          const assemble::AANeighborListContainer current_chain_neighbor_list_container
          (
            amino_acids_current_chain,
            ols.GetDistanceCutoff(),
            ols.GetMinimalSequenceSeparation(),
            false
          );

          const assemble::AANeighborListContainer &desired_neighbor_list_container
          (
            m_ChainOption == e_OneChain
            ? current_chain_neighbor_list_container
            : all_chain_neighbor_list_container
          );

          // iterate through current chain
          for
          (
            util::SiPtrVector< const biol::AABase>::const_iterator
              aa_itr( amino_acids_current_chain.Begin()), aa_itr_end( amino_acids_current_chain.End());
            aa_itr != aa_itr_end;
            ++aa_itr
          )
          {
            // proceed only if there are coordinates given and the aa is a natural aa
            if
            (
              !( *aa_itr)->GetFirstSidechainAtom().GetCoordinates().IsDefined() ||
              !( *aa_itr)->GetType()->IsNaturalAminoAcid()
            )
            {
              continue;
            }

            // get current environment type
            const biol::EnvironmentType current_environment_type
            (
              sp_membrane.IsDefined() && m_SplitEnvironment
              ? sp_membrane->DetermineEnvironmentType( ( *aa_itr)->GetFirstSidechainAtom().GetCoordinates())
              : biol::GetEnvironmentTypes().e_Solution
            );

            // calculate desired ols
            const double &desired_ols
            (
              ols( desired_neighbor_list_container.Find( ( **aa_itr))->second)
            );

            desired_histograms( current_environment_type)( ( *aa_itr)->GetType()).PushBack( desired_ols);

          } // end of iterating through all aas in current chain
        } // end of iterating through all chains
      } // end of iterating through all protein models

      // write output
      std::stringstream stream;

      stream << ols.GetThresholdRange() << '\n';
      stream << ols.GetMinimalSequenceSeparation() << '\n';

      // iterate over environment types
      for
      (
        biol::EnvironmentTypes::const_iterator
          ent_itr( biol::GetEnvironmentTypes().Begin()), ent_itr_end( biol::GetEnvironmentTypes().End());
        ent_itr != ent_itr_end;
        ++ent_itr
      )
      {
        // only interested in soluble environment type
        if( !m_SplitEnvironment)
        {
          if( ( *ent_itr) != biol::GetEnvironmentTypes().e_Solution)
          {
            continue;
          }
        }
        else
        {
          stream << *ent_itr << '\n';
        }

        // for current environment type iterate over all amino acid types
        for
        (
          biol::AATypes::const_iterator
            aa_itr( biol::GetAATypes().Begin()),
            aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          stream << ( *aa_itr)->GetOneLetterCode() << '\n';
          stream << desired_histograms( *ent_itr)( *aa_itr);
        }
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer OLS::GetSerializer() const
    {

      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes neighbor count statistics."
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "sequence_exclusion",
        "minimal sequence separation to be considered in the calculation",
        io::Serialization::GetAgent( &m_SequenceExclusion),
        "2"
      );

      parameters.AddInitializer
      (
        "radius",
        "radius for the overlapping sphere algorithm",
        io::Serialization::GetAgent( &m_Radius),
        "4.75"
      );

      parameters.AddInitializer
      (
        "chain_option",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_ChainOption),
        GetChainOptionName( e_OneChain)
      );
      return parameters;
    }
  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_phipsi.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram_2d.h"
#include "sspred/bcl_sspred_ci_phi_psi.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PhiPsi::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new PhiPsi())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &PhiPsi::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Histogram2D",
        GetStaticClassName< PhiPsi::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &PhiPsi::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "phi_psi_angles_by_sstype.histogram2D",
        GetStaticClassName< PhiPsi::OutputOption>()
      };
      return s_output_filename_extensions[ size_t(OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PhiPsi::PhiPsi() :
      m_OutputOption( e_Histogram2D),
      m_NumberOfBins( 12),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    PhiPsi *PhiPsi::Clone() const
    {
      return new PhiPsi( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &PhiPsi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &PhiPsi::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns the number of bins of Phi/Psi dihedral historgams
    //! @return the number of bins of Phi/Psi dihedral histograms
    const size_t &PhiPsi::GetNumberOfBins() const
    {
      return m_NumberOfBins;
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &PhiPsi::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &PhiPsi::GetAlias() const
    {
      static std::string s_name( "PhiPsi");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string PhiPsi::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // statistics of phi/psi dihedrals
      storage::Vector< storage::Vector< storage::Vector< math::Histogram2D> > > phi_psi_sstypes_histograms
      (
        2,
        storage::Vector< storage::Vector< math::Histogram2D> >
        (
          biol::GetSSTypes().GetEnumCount(),
          storage::Vector< math::Histogram2D>
          (
            biol::GetAATypes().GetEnumCount(),
            math::Histogram2D
            (
              storage::VectorND< 2, double>( -math::g_Pi, -math::g_Pi),
              storage::VectorND< 2, double>( 2 * math::g_Pi / m_NumberOfBins, 2 * math::g_Pi / m_NumberOfBins),
              storage::VectorND< 2, size_t>( m_NumberOfBins, m_NumberOfBins)
            )
          )
        )
      );

      // iterate over the protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        util::ShPtr< assemble::ProteinModel> sp_model( *protein_model_itr);
        sspred::CIPhiPsi().Calculate( *sp_model, true);

        // get the current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get membrane for current protein model
        const util::ShPtr< biol::Membrane> sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        // get all chains in current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // iterate over all chains in current protein model
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator
            chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            BCL_MessageStd( "Skip undesired chain: " + util::Format()( ( *chain_itr)->GetChainID()));
            continue;
          }

          // get all sses in current chain
          const util::SiPtrVector< const assemble::SSE> &all_sses( ( *chain_itr)->GetSSEs());

          // iterate over all sse in current chain
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            // skip undefined sses
            if( !( *sse_itr)->GetType().IsDefined())
            {
              BCL_MessageStd( "Skip undefined sse: " + ( *sse_itr)->GetType().GetName());
              continue;
            }

            // get all amino acids in current sse
            const util::ShPtrVector< biol::AABase> &all_amino_acids( ( *sse_itr)->GetData());

            // iterate over all amino acid in current sse
            for
            (
              util::ShPtrVector< biol::AABase>::const_iterator
                aa_itr( all_amino_acids.Begin()), aa_itr_end( all_amino_acids.End());
              aa_itr != aa_itr_end;
              ++aa_itr
            )
            {
              // skip unnatural amino acids or those that are first or last in the sse
              if
              (
                !( *aa_itr)->GetType()->IsNaturalAminoAcid() ||
                aa_itr == all_amino_acids.Begin() ||
                aa_itr == aa_itr_end - 1
              )
              {
                continue;
              }

              // calculate phi/psi using C atom from preceding residue and N atom from following residue
              storage::VectorND< 2, double> phi_psi_vector
              (
                ( *aa_itr)->CalculatePhi( ( *( aa_itr - 1))->GetAtom( biol::GetAtomTypes().C)),
                ( *aa_itr)->CalculatePsi( ( *( aa_itr + 1))->GetAtom( biol::GetAtomTypes().N))
              );

              if( m_OutputOption == e_Histogram2D)
              {
                auto sspre( ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_CIPhiPsi));
                biol::EnvironmentType env_type
                (
                  sp_membrane.IsDefined() && sspre.IsDefined()
                  ? sspre->GetOneStateTMPrediction()->GetReducedType()
                  : biol::GetEnvironmentTypes().e_Solution
                );
                phi_psi_sstypes_histograms( env_type == biol::GetEnvironmentTypes().e_MembraneCore ? 1 : 0)
                ( ( *sse_itr)->GetType())
                ( ( *aa_itr)->GetType()).PushBack( phi_psi_vector);
              }
            } // end of iteration over amino acid residues
          } // end of iteration over sses
        } // end of iteration over chains
      } // end of iteration over protein models

      // write statistics
      std::stringstream stream;
      if( m_OutputOption == e_Histogram2D)
      {
        // iterate over sse types
        for( size_t in_membrane( 0), im_max( 2); in_membrane < im_max; ++in_membrane)
        {
          stream << ( in_membrane ? "MEMBRANE" : "SOLUTION") << '\n';
          for
          (
            biol::SSTypes::const_iterator
              sse_type_itr( biol::GetSSTypes().Begin()), sse_type_itr_end( biol::GetSSTypes().COIL.GetIterator() + 1);
            sse_type_itr != sse_type_itr_end;
            ++sse_type_itr
          )
          {
            // write sse type
            stream << *sse_type_itr << '\n';

            // iterate over aa types
            for
            (
              biol::AATypes::const_iterator
                aa_itr( biol::GetAATypes().Begin()),
                aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
              aa_itr != aa_itr_end;
              ++aa_itr
            )
            {
              // write one letter code of aa
              stream << ( *aa_itr)->GetOneLetterCode() << '\n';
              // write histograms
              stream << phi_psi_sstypes_histograms( in_membrane)( *sse_type_itr)( *aa_itr);
            }
          }
        }
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PhiPsi::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes phi/psi dihedral statistics."
      );

      parameters.AddInitializer
      (
        "resolution",
        "number of bins of the phi/psi histogram",
        io::Serialization::GetAgent( &m_NumberOfBins),
        "12"
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "Histogram2D"
      );

      return parameters;
    }

  } // namespace scorestat
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
// include header for this class
#include "scorestat/bcl_scorestat_protein_model_packing.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "biol/bcl_biol_aa_back_bone.h"
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
#include "score/bcl_score_aa_pair_hi_res_clash.h"
#include "score/bcl_score_aa_pair_sidechain_interaction.h"
#include "storage/bcl_storage_table.h"
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
    const util::SiPtr< const util::ObjectInterface> ProteinModelPacking::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new ProteinModelPacking())
    );

    //! @brief Category as string
    //! @param CATEGORY the desired category
    //! @return the string for the category
    const std::string &ProteinModelPacking::GetCategoryName( const Category &CATEGORY)
    {
      static std::string s_names[ s_NumberCategories + 1] =
      {
        "AdjacentInContact",
        "AdjacentNotInContact",
        "AdjacentParallel",
        "AdjacentAntiParallel",
        "OneSSEApartParallel",
        "OneSSEApartAntiParallel",
        "TwoSSEApartParallel",
        "TwoSSEApartAntiParallel",
        "FarApartParallel",
        "FarApartAntiParallel",
        "WeakInteraction",
        "ModerateInteraction",
        "StrongInteraction",
        "BackgroundWeakInteraction",
        "BackgroundModerateInteraction",
        "BackgroundStrongInteraction",
        "BackgroundParallel",
        "BackgroundAntiParallel",
        GetStaticClassName< Category>()
      };
      return s_names[ CATEGORY];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelPacking::ProteinModelPacking
    (
      const double &MIN_INTERACTION_DISTANCE,
      const size_t &MIN_ATOMS_IN_CONTACT
    ) :
      m_InteractionDistance( MIN_INTERACTION_DISTANCE),
      m_MinAtomsInContact( MIN_ATOMS_IN_CONTACT)
    {
    }

    //! @brief virtual copy constructor
    ProteinModelPacking *ProteinModelPacking::Clone() const
    {
      return new ProteinModelPacking( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelPacking::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &ProteinModelPacking::GetOutFilePostfix() const
    {
      static const std::string s_name( "sse_packing_type");
      return s_name;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ProteinModelPacking::GetAlias() const
    {
      static const std::string s_name( "SSEPackingType");
      return s_name;
    }

    //! @brief hash a given packing. The given string is intended to be fast to hash, not necessarily easily readable
    void ProteinModelPacking::AddPackingType
    (
      const assemble::SSEGeometryPacking &PACKING,
      const bool &IS_IN_CONTACT,
      const size_t &N_SSES_APART,
      const bool &IS_BACKGROUND,
      const bool &ORIENTATION_COULD_BE_OPPOSITE,
      linal::VectorInterface< size_t> &COUNTS
    )
    {
      if( !IS_BACKGROUND)
      {
        // adjacency contacts. Bins 0 - 1
        // Orientation contacts. Bins 2 - 5 (Adjacent Parallel, Adjacent Anti-Parallel, Non-Adj Parallel, Non-Adj Anti-Parallel)
        if( !N_SSES_APART)
        {
          ++COUNTS( IS_IN_CONTACT ? e_AdjacentInContact : e_AdjacentNotInContact);
        }
        if( !IS_IN_CONTACT)
        {
          return;
        }
        if( ORIENTATION_COULD_BE_OPPOSITE)
        {
          ++COUNTS
          (
            std::min( N_SSES_APART, size_t( 3)) * 2 + size_t( e_AdjacentParallel)
            + ( PACKING.GetOrientation() == assemble::SSEGeometryPacking::e_AntiParallel ? 1 : 0)
          );
        }
        ++COUNTS
        (
          PACKING.GetInteractionWeight() < 0.5
          ? e_WeakInteraction
          : (
              PACKING.GetInteractionWeight() < 0.95
              ? e_ModerateInteraction
              : e_StrongInteraction
            )
        );
      }
      else
      {
        ++COUNTS
        (
          PACKING.GetInteractionWeight() < 0.5
          ? e_BackgroundWeakInteraction
          : (
              PACKING.GetInteractionWeight() < 0.95
              ? e_BackgroundModerateInteraction
              : e_BackgroundStrongInteraction
            )
        );
        ++COUNTS
        (
          PACKING.GetOrientation() == assemble::SSEGeometryPacking::e_Parallel
          ? e_BackgroundParallel
          : e_BackgroundAntiParallel
        );
      }
    }

  //////////////
  // operator //
  //////////////

    //! @brief add statistics about the given protein and chains to the internal table and histograms
    //! @param ENSEMBLE the protein ensemble the statistics of which to be computed
    std::string ProteinModelPacking::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // min probability to score two residues as being an almost-certain contact
      const double min_p( 0.04);

      // initializes maps for storing histograms of distances between all types of amino acid pairs
      storage::Vector< linal::Vector< size_t> > packing_info_by_contact_type
      (
        contact::GetTypes().GetEnumCount(),
        linal::Vector< size_t>( size_t( s_NumberCategories), size_t( 0))
      );
      storage::Vector< linal::Vector< size_t> > contact_type_by_ss_types
      (
        size_t( 4),
        linal::Vector< size_t>( size_t( contact::GetTypes().GetEnumCount()), size_t( 0))
      );
      storage::Vector< linal::Vector< size_t> > contact_type_by_ss_types_bg
      (
        size_t( 4),
        linal::Vector< size_t>( size_t( contact::GetTypes().GetEnumCount()), size_t( 0))
      );
      const long max_helix_helix_contacts( 8), max_strand_strand_contacts( 6);
      const long max_helix_strand_contacts( 6), max_strand_helix_contacts( 3);

      linal::Matrix< double> helix_adjacency_counts( max_helix_helix_contacts + 1, max_helix_strand_contacts + 1, double( 0.0));
      linal::Matrix< double> strand_adjacency_counts( max_strand_helix_contacts + 1, max_strand_strand_contacts + 1, double( 0.0));
      linal::Matrix< double> helix_adjacency_plaus( helix_adjacency_counts), strand_adjacency_plaus( strand_adjacency_counts);
      storage::Vector< double> counts_contacts_possible( size_t( 4));
      storage::Vector< double> counts_contacts_actual( size_t( 4));
      // for each sse, compute average contact type of contacting SSEs, weighted by # of AAs in contact
      // Also compute P(Helix|Helix) In contact, P(Strand|Strand)
      // # of helix-strand contacts / (min( # of helix * min(# of strands,4), min(# of helices,2)*# of strands))
      // # of helix-helix contacts * 2 / (# of helices * min(# of helices - 1,4))
      // # of strand-strand contacts * 2 / (# of strands * min(# of strand - 1,4))

      storage::Set< biol::AtomType> types_of_interest( biol::GetAtomTypes().GetFirstSidechainAtomTypes());

      const std::string sse_name[ 2] = { "HELIX", "STRAND"};

      size_t n_runs( 10000);
      score::AAPairHiResClash clash;
      for( int config( 0); config < 4; ++config)
      {
        const size_t hs_a( config & 1 ? 1 : 0);
        const size_t hs_b( config & 2 ? 1 : 0);
        BCL_MessageStd( "Computing background for " + sse_name[ hs_a] + " and " + sse_name[ hs_b]);

        for( size_t run_n( 0), n_runs_l( n_runs); run_n < n_runs_l; ++run_n)
        {
          const size_t strand_l( 6), helix_l( 10);
          storage::VectorND< 2, storage::VectorND< 3, assemble::SSE> > helix_strand_res;
          int o_id( 1);
          int pdb_atom_id( 1);
          for( size_t o( 0); o < size_t( 3); ++o, o_id += 2)
          {
            util::ShPtrVector< biol::AABase> strand_aas( strand_l);
            for( size_t i( 0); i < strand_l; ++i, ++o_id)
            {
              double ran_choice( random::GetGlobalRandom().Random( 1.0));
              biol::AAType type( biol::GetAATypes().ALA);
              ran_choice -= type->GetAAProperty( biol::AATypeData::e_NaturalPrevalence);
              for( int type_x( 1); ran_choice > 0.0; type_x += 1)
              {
                type = biol::AAType( type_x);
                ran_choice -= type->GetAAProperty( biol::AATypeData::e_NaturalPrevalence);
              }
              util::ShPtr< biol::AAData> aadata( new biol::AAData( type, o_id, o_id, 'A'));
              strand_aas( i) = util::ShPtr< biol::AABase>( new biol::AABackBone( aadata));
              strand_aas( i)->SetAtom( biol::Atom( biol::GetAtomTypes().N, pdb_atom_id++));
              strand_aas( i)->SetAtom( biol::Atom( biol::GetAtomTypes().C, pdb_atom_id++));
              strand_aas( i)->SetAtom( biol::Atom( biol::GetAtomTypes().CA, pdb_atom_id++));
              strand_aas( i)->SetAtom( biol::Atom( biol::GetAtomTypes().O, pdb_atom_id++));
              strand_aas( i)->SetAtom( biol::Atom( type->GetFirstSidechainAtomType(), pdb_atom_id++));
            }
            o_id += 4;
            util::ShPtrVector< biol::AABase> helix_aas( helix_l);
            for( size_t j( 0); j < helix_l; ++j, ++o_id)
            {
              double ran_choice( random::GetGlobalRandom().Random( 1.0));
              biol::AAType type( biol::GetAATypes().ALA);
              ran_choice -= type->GetAAProperty( biol::AATypeData::e_NaturalPrevalence);
              for( int type_x( 1); ran_choice > 0.0; type_x += 1)
              {
                type = biol::AAType( type_x);
                ran_choice -= type->GetAAProperty( biol::AATypeData::e_NaturalPrevalence);
              }
              util::ShPtr< biol::AAData> aadata( new biol::AAData( type, o_id, o_id, 'A'));
              helix_aas( j) = util::ShPtr< biol::AABase>( new biol::AABackBone( aadata));
              helix_aas( j)->SetAtom( biol::Atom( biol::GetAtomTypes().N, pdb_atom_id++));
              helix_aas( j)->SetAtom( biol::Atom( biol::GetAtomTypes().C, pdb_atom_id++));
              helix_aas( j)->SetAtom( biol::Atom( biol::GetAtomTypes().CA, pdb_atom_id++));
              helix_aas( j)->SetAtom( biol::Atom( biol::GetAtomTypes().O, pdb_atom_id++));
              helix_aas( j)->SetAtom( biol::Atom( type->GetFirstSidechainAtomType(), pdb_atom_id++));
            }
            o_id += 4;
            helix_strand_res( 0)( o) = assemble::SSE( biol::AASequence( helix_aas), biol::GetSSTypes().HELIX);
            helix_strand_res( 1)( o) = assemble::SSE( biol::AASequence( strand_aas), biol::GetSSTypes().STRAND);
            helix_strand_res( 0)( o).SetToIdealConformationAtOrigin();
            helix_strand_res( 1)( o).SetToIdealConformationAtOrigin();
          }
          assemble::SSE &sse_a( helix_strand_res( hs_a)( 0));
          assemble::SSE &sse_b( helix_strand_res( hs_b)( 1));

          sse_a.SetToIdealConformationAtOrigin();
          sse_b.SetToIdealConformationAtOrigin();
          sse_a.Rotate( math::RotationMatrix3D().SetRand());
          sse_b.Rotate( math::RotationMatrix3D().SetRand());
          sse_a.Translate( linal::Vector3D().SetRandomTranslation( 20.0));
          sse_b.Translate( linal::Vector3D().SetRandomTranslation( 20.0));
          if( clash( sse_a, sse_b))
          {
            --run_n;
            continue;
          }

          assemble::VoxelGridAA interactions_detector( m_InteractionDistance);
          util::SiPtrVector< const assemble::SSE> structured_sses( util::SiPtrVector< const assemble::SSE>::Create( sse_a, sse_b));
          util::SiPtrVector< const biol::AABase> aas( sse_a.Begin(), sse_a.End());
          aas.Append( util::SiPtrVector< const biol::AABase>( sse_b.Begin(), sse_b.End()));
          linal::Matrix< float> interactions_matrix
          (
            interactions_detector.GetSSEInteractionMatrix( structured_sses, aas, 2, m_InteractionDistance, false, min_p, false)
          );

          if( interactions_matrix( 0, 1) < m_MinAtomsInContact)
          {
            --run_n;
            continue;
          }
          // at this point we know that we have a Beta-Alpha-Beta with no intervening SSEs. Next, check that the strands
          // point in the same general direction
          assemble::SSEGeometryPacking packing_b( sse_a, sse_b, 0.0);
          if( !packing_b.GetContactType().IsDefined())
          {
            continue;
          }

          AddPackingType( packing_b, true, 1, true, true, packing_info_by_contact_type( packing_b.GetContactType()));
          AddPackingType( packing_b, true, 1, true, true, packing_info_by_contact_type( contact::GetTypes().Reverse( packing_b.GetContactType())));
          ++contact_type_by_ss_types_bg( config)( packing_b.GetContactType());
          if( !( run_n % ( n_runs_l / 100)))
          {
            util::GetLogger().LogStatus( util::Format()( run_n * 100 / n_runs_l) + "% complete");
          }
        }
      }

      // Create a matrix that will hold the x,y, and z coordinates for terminii of the first & second strands
      for( auto itr_ensemble( ENSEMBLE.Begin()), itr_ensemble_end( ENSEMBLE.End()); itr_ensemble != itr_ensemble_end; ++itr_ensemble)
      {
        // get current protein model
        util::ShPtr< assemble::ProteinModel> sp_current_model( *itr_ensemble);
        const assemble::ProteinModel &protein_model( *sp_current_model);

        // get current pdb name before all the loops start
        const util::ShPtr< util::Wrapper< std::string> > &model_name_ptr
        (
          protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
        );
        const std::string model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");

        // instantiate sequences CaCb from pdb for each chain one
        BCL_MessageStd( "computing packing for " + model_name);

        for
        (
          auto itr_chain( ( *itr_ensemble)->GetChains().Begin()), itr_chain_end( ( *itr_ensemble)->GetChains().End());
          itr_chain != itr_chain_end;
          ++itr_chain
        )
        {
          // need at least three sses in the chain
          if( ( *itr_chain)->GetNumberSSEs() < 2)
          {
            continue;
          }
          // collect all non-coil sse's
          const util::SiPtrVector< const assemble::SSE> structured_sses
          (
            ( *itr_chain)->GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
          );
          const long n_helices( ( *itr_chain)->GetSSEs( biol::GetSSTypes().HELIX).GetSize());
          const long n_strands( ( *itr_chain)->GetSSEs( biol::GetSSTypes().STRAND).GetSize());

          const double helix_exp
          (
            double( n_helices)
            / std::max
              (
                1.0,
                double
                (
                  std::min( n_helices, long( max_helix_helix_contacts + 1)) *
                  std::min( n_strands + 1, long( max_helix_strand_contacts + 1))
                )
              )
          );
          const double strand_exp
          (
            double( n_strands)
            / std::max
              (
                1.0,
                double
                (
                  std::min( n_strands, long( max_strand_strand_contacts + 1)) *
                  std::min( n_helices + 1, long( max_strand_helix_contacts + 1))
                )
              )
          );

          for( long i( 0); i <= n_helices; ++i)
          {
            for( long j( 0); j <= n_strands; ++j)
            {
              if( i < n_helices && i <= max_helix_helix_contacts && j <= max_helix_strand_contacts)
              {
                helix_adjacency_plaus( i, j) += helix_exp;
              }
              if( j < n_strands && i <= max_strand_helix_contacts && j <= max_strand_strand_contacts)
              {
                strand_adjacency_plaus( i, j) += strand_exp;
              }
            }
          }

          counts_contacts_possible( 0) += n_helices * double( std::min( long( 6), std::max( long( n_helices - 1), long( 0))));
          counts_contacts_possible( 3) += n_strands * double( std::min( long( 4), std::max( long( n_strands - 1), long( 0))));
          const long expected_n_helix_strand
          (
            std::max( double( n_helices * std::min( n_strands, long( 4))), double( std::min( n_helices, long( 4)) * n_strands))
          );
          counts_contacts_possible( 1) += expected_n_helix_strand;
          counts_contacts_possible( 2) += expected_n_helix_strand;

          if( structured_sses.GetSize() < 2)
          {
            continue;
          }

          util::SiPtrVector< const biol::Atom> si_atoms;

          for( auto itr_sse( structured_sses.Begin()), itr_sse_end( structured_sses.End()); itr_sse != itr_sse_end; ++itr_sse)
          {
            si_atoms.Append( ( *itr_sse)->GetAtoms( types_of_interest));
          }
          assemble::VoxelGridAA interactions_detector( m_InteractionDistance);
          linal::Matrix< float> interactions_matrix
          (
            interactions_detector.GetSSEInteractionMatrix
            (
              structured_sses,
              ( *itr_chain)->GetAminoAcids(),
              2,
              m_InteractionDistance,
              false,
              min_p,
              false
            )
          );
          for( size_t i( 0), nr_sse( interactions_matrix.GetNumberCols()); i < nr_sse; ++i)
          {
            size_t strand_count( 0), helix_count( 0);
            for( size_t j( 0); j < nr_sse; ++j)
            {
              if( interactions_matrix( i, j) >= m_MinAtomsInContact)
              {
                ++( structured_sses( j)->GetType() == biol::GetSSTypes().STRAND ? strand_count : helix_count);
              }
            }
            linal::Matrix< double> &ref_mat( structured_sses( i)->GetType() == biol::GetSSTypes().STRAND ? strand_adjacency_counts : helix_adjacency_counts);
            if( helix_count < ref_mat.GetNumberRows() && strand_count < ref_mat.GetNumberCols())
            {
              ref_mat( helix_count, strand_count) += 1.0;
            }
          }

          // handle SSE-pair interaction-weight and parallel/anti-parallel bias
          size_t a( 0);
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              itr_a( structured_sses.Begin()),
              itr_end( structured_sses.End());
            itr_a != itr_end;
            ++itr_a, ++a
          )
          {
            size_t b( a + 1);
            const size_t a_type( ( *itr_a)->GetType() == biol::GetSSTypes().HELIX ? 0 : 1);
            for
            (
              util::SiPtrVector< const assemble::SSE>::const_iterator itr_b( itr_a + 1);
              itr_b != itr_end;
              ++itr_b, ++b
            )
            {
              const size_t b_type( ( *itr_b)->GetType() == biol::GetSSTypes().HELIX ? 0 : 2);
              const bool is_in_contact( interactions_matrix( a, b) >= m_MinAtomsInContact);
              const assemble::SSEGeometryPacking::Orientation ori
              (
                assemble::SSEGeometryPacking::OrientationFromSSEs( ( **itr_a), ( **itr_b))
              );
              // for every pair of SSEs
              assemble::SSEGeometryPacking packing_ab( **itr_a, **itr_b, 0.0);
              if( !packing_ab.GetContactType().IsDefined())
              {
                continue;
              }
              const contact::Type reverse_type( contact::GetTypes().Reverse( packing_ab.GetContactType()));
              if( !is_in_contact && itr_b != itr_a + 1)
              {
                continue;
              }

              const double max_loop_length
              (
                ( ( *itr_b)->GetFirstAA()->GetPdbID() - ( *itr_a)->GetLastAA()->GetPdbID()) * 2.56 + 2.11
              );
              const bool could_be_opposite_ori
              (
                ( *itr_b)->GetChainID() != ( *itr_a)->GetChainID()
                ||
                std::min
                (
                  linal::Distance( ( *itr_b)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_a)->GetLastAA()->GetCA().GetCoordinates()),
                  std::min
                  (
                    linal::Distance( ( *itr_b)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_b)->GetFirstAA()->GetCA().GetCoordinates()),
                    linal::Distance( ( *itr_a)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_a)->GetFirstAA()->GetCA().GetCoordinates())
                  ) + 2.34
                ) < max_loop_length
              );
//              if( is_in_contact)
//              {
//                BCL_MessageStd
//                (
//                  std::string( itr_b == itr_a + 1 ? "Adj " : "NonAdj ")
//                  + std::string( could_be_opposite_ori ? " CouldSwitch " : " Constrained ")
//                  + std::string( is_parallel ? " Parallel " : " Antiparallel ")
//                  + util::Format()( max_loop_length) + " "
//                  + util::Format()( linal::Distance( ( *itr_b)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_a)->GetLastAA()->GetCA().GetCoordinates()))
//                  + " " + util::Format()( linal::Distance( ( *itr_b)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_b)->GetFirstAA()->GetCA().GetCoordinates()) + 2.34)
//                  + " " + util::Format()( linal::Distance( ( *itr_a)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_a)->GetFirstAA()->GetCA().GetCoordinates()) + 2.34)
//                  + " " + ( *itr_a)->GetIdentification() + " " + ( *itr_b)->GetIdentification()
//                );
//              }
              const size_t sse_separation( itr_b - itr_a - 1);
              AddPackingType
              (
                packing_ab,
                is_in_contact,
                sse_separation,
                false,
                could_be_opposite_ori,
                packing_info_by_contact_type( packing_ab.GetContactType())
              );
              AddPackingType
              (
                packing_ab,
                is_in_contact,
                sse_separation,
                false,
                could_be_opposite_ori,
                packing_info_by_contact_type( reverse_type)
              );
              if( is_in_contact)
              {
                counts_contacts_actual( b_type | a_type) += 1;
                if( a_type != b_type)
                {
                  counts_contacts_actual( ( a_type ? 2 : 0) + ( b_type ? 1 : 0)) += 1;
                }
              }
              ++contact_type_by_ss_types( a_type + b_type)( packing_ab.GetContactType());
              ++contact_type_by_ss_types( ( a_type ? 2 : 0) + ( b_type ? 1 : 0))( reverse_type);
            }
          }
        }
      }

      // write statistics
      std::ostringstream stream;

      stream << this->GetString() << '\n';

      const size_t n_contact_types( contact_type_by_ss_types( 0).GetSize());
      stream << "Type\t";
      for
      (
        auto itr_ct( contact::GetTypes().Begin()), itr_ct_end( contact::GetTypes().End());
        itr_ct != itr_ct_end;
        ++itr_ct
      )
      {
        stream << itr_ct->GetName() << '\t';
      }

      stream << '\n';
      stream << "LogOdds\t";
      for( size_t ct( 0); ct < n_contact_types; ++ct)
      {
        int config( 0);
        while( config < 4 && !contact_type_by_ss_types_bg( config)( ct))
        {
          ++config;
        }
        const linal::Vector< size_t> &counts_real( contact_type_by_ss_types( config));
        const linal::Vector< size_t> &counts_bg( contact_type_by_ss_types_bg( config));
        const size_t sum_real( counts_real.Sum());
        const size_t sum_bg( counts_bg.Sum());
        const double counts_ratio( double( sum_real) / double( sum_bg));
        stream << -std::log( ( double( counts_real( ct)) / double( counts_bg( ct))) / counts_ratio) << '\t';
      }
      stream << '\n';
      stream << "\nContactType\t";
      for( size_t i( 0); i < s_NumberNaturalCategories; ++i)
      {
        stream << CategoryEnum( Category( i)).GetString() << '\t';
      }
      stream << '\n';
      for( size_t ct( 0); ct < n_contact_types; ++ct)
      {
        stream << contact::Type( ct).GetName() << '\t';
        const linal::Vector< size_t> &contact_counts( packing_info_by_contact_type( ct));
        const double adjacent_contact( contact_counts( e_AdjacentParallel) + contact_counts( e_AdjacentAntiParallel));
        const double one_sse_contact( contact_counts( e_OneSSEApartParallel) + contact_counts( e_OneSSEApartAntiParallel));
        const double two_sse_contact( contact_counts( e_TwoSSEApartParallel) + contact_counts( e_TwoSSEApartAntiParallel));
        const double three_sse_contact( contact_counts( e_ThreeOrMoreSSEApartParallel) + contact_counts( e_ThreeOrMoreSSEApartAntiParallel));
        const double real_total
        (
          contact_counts( e_WeakInteraction)
          + contact_counts( e_ModerateInteraction)
          + contact_counts( e_StrongInteraction)
        );
        const double bg_total
        (
          contact_counts( e_BackgroundWeakInteraction)
          + contact_counts( e_BackgroundModerateInteraction)
          + contact_counts( e_BackgroundStrongInteraction)
        );

        const double bg_parallel_propensity
        (
          double( contact_counts( e_BackgroundParallel)) / bg_total
        );
        const double bg_antiparallel_propensity
        (
          double( contact_counts( e_BackgroundAntiParallel)) / bg_total
        );
        const double real_bg_ratio( real_total / bg_total);
        const double adjacent_total( contact_counts( e_AdjacentInContact) + contact_counts( e_AdjacentNotInContact));
        stream << -std::log( double( contact_counts( e_AdjacentInContact) * 2) / adjacent_total) << '\t'
               << -std::log( double( contact_counts( e_AdjacentNotInContact) * 2) / adjacent_total) << '\t'
               << -std::log( double( contact_counts( e_AdjacentParallel)) / adjacent_contact / bg_parallel_propensity) << '\t'
               << -std::log( double( contact_counts( e_AdjacentAntiParallel)) / adjacent_contact / bg_antiparallel_propensity) << '\t'
               << -std::log( double( contact_counts( e_OneSSEApartParallel) / one_sse_contact / bg_parallel_propensity)) << '\t'
               << -std::log( double( contact_counts( e_OneSSEApartAntiParallel) / one_sse_contact / bg_antiparallel_propensity)) << '\t'
               << -std::log( double( contact_counts( e_TwoSSEApartParallel) / two_sse_contact / bg_parallel_propensity)) << '\t'
               << -std::log( double( contact_counts( e_TwoSSEApartAntiParallel) / two_sse_contact / bg_antiparallel_propensity)) << '\t'
               << -std::log( double( contact_counts( e_ThreeOrMoreSSEApartParallel) / three_sse_contact / bg_parallel_propensity)) << '\t'
               << -std::log( double( contact_counts( e_ThreeOrMoreSSEApartAntiParallel) / three_sse_contact / bg_antiparallel_propensity)) << '\t'
               << -std::log( double( contact_counts( e_WeakInteraction) + real_bg_ratio) / double( contact_counts( e_BackgroundWeakInteraction) + 1) / real_bg_ratio) << '\t'
               << -std::log( double( contact_counts( e_ModerateInteraction) + real_bg_ratio) / double( contact_counts( e_BackgroundModerateInteraction) + 1) / real_bg_ratio) << '\t'
               << -std::log( double( contact_counts( e_StrongInteraction) + real_bg_ratio) / double( contact_counts( e_BackgroundStrongInteraction) + 1) / real_bg_ratio) << '\n';
      }
      stream << '\n';
      stream << "\n\n# The following information is not currently use but is collected to aid in the development of further scores\n"
             << "SS-PairType\tMaxPairingsPossible\tActualNumber\n"
             << "HelixHelix\t" << counts_contacts_possible( 0) << '\t' << counts_contacts_actual( 0) << '\n'
             << "HelixStrand\t" << counts_contacts_possible( 1) << '\t' << counts_contacts_actual( 1) << '\n'
             << "StrandStrand\t" << counts_contacts_possible( 3) << '\t' << counts_contacts_actual( 3) << '\n';

      stream << "\n# Adjacency Counts for strands. Rows - # of helices in contact. Columns - # of strands in contact\n";
      stream << "\n" << strand_adjacency_counts;
      stream << "\nAdjacency Counts for helices. Rows - # of helices in contact. Columns - # of strands in contact\n"
             << "\n" << helix_adjacency_counts;
      auto itr_plaus( strand_adjacency_plaus.Begin());
      for( auto itr( strand_adjacency_counts.Begin()), itr_end( strand_adjacency_counts.End()); itr != itr_end; ++itr, ++itr_plaus)
      {
        *itr /= *itr_plaus;
      }
      itr_plaus = helix_adjacency_plaus.Begin();
      for( auto itr( helix_adjacency_counts.Begin()), itr_end( helix_adjacency_counts.End()); itr != itr_end; ++itr, ++itr_plaus)
      {
        *itr /= *itr_plaus;
      }
      stream << "\n#Adjacency Propensity for strands (# times / # times expected from background w/ same helix-strand distribution) for helices.\n"
             << "\n#Rows - # of helices in contact. Columns - # of strands in contact\n";
      stream << "\n" << strand_adjacency_counts << "\n";
      stream << "\n#Adjacency propensity for helices (# times / # times expected from background w/ same helix-strand distribution) for helices.\n"
              << "\n#Rows - # of helices in contact. Columns - # of strands in contact\n";
      stream << helix_adjacency_counts;

      return stream.str();
    } // end of operator

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelPacking::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Computes packing type, distance, interaction weight entropies");

      parameters.AddInitializer
      (
        "distance",
        "max distance between atoms to be considered in contact",
        io::Serialization::GetAgent( &m_InteractionDistance),
        util::Format()( ProteinModelPacking().m_InteractionDistance)
      );
      parameters.AddInitializer
      (
        "min sc atoms",
        "minimum side chains atoms betwen contacting SSEs",
        io::Serialization::GetAgent( &m_MinAtomsInContact),
        util::Format()( ProteinModelPacking().m_MinAtomsInContact)
      );
      return parameters;
    }

    //! @brief read in the SSPair to contact entropy table into a vector with
    //! @param STREAM Input stream to read from
    linal::Vector< double> ProteinModelPacking::ReadSSPairToContactEntropies( std::istream &STREAM) const
    {
      storage::Table< double> table;
      table.ReadFormatted( STREAM);
      BCL_Assert( table.GetNumberRows() == size_t( 1), "Table should have size 4!");
      linal::Vector< double> contact_type_entropies( contact::GetTypes().GetEnumCount());
      const size_t n_cols( table.GetNumberColumns());
      auto itr_col( table.Begin()->Second().GetData().Begin());
      for( size_t col_number( 0); col_number < n_cols; ++col_number, ++itr_col)
      {
        contact_type_entropies( contact::Type( table.GetHeader()( col_number)).GetIndex()) = *itr_col;
      }
      return contact_type_entropies;
    }

    //! @brief read in the entropy table into a vec of vecs. Rows are indexed by contact type, columns by Categories Enum
    //! @param STREAM Input stream to read from
    storage::Vector< linal::Vector< double> > ProteinModelPacking::ReadContactTypeEntropies( std::istream &STREAM) const
    {
      storage::Table< double> table;
      table.ReadFormatted( STREAM);
      storage::Vector< linal::Vector< double> > contact_entropies
      (
        contact::GetTypes().GetEnumCount(),
        linal::Vector< double>( size_t( s_NumberNaturalCategories), double( 0.0))
      );
      const size_t n_cols( table.GetNumberColumns());
      for( auto itr( table.Begin()), itr_end( table.End()); itr != itr_end; ++itr)
      {
        const size_t row_number( contact::Type( itr->First()).GetIndex());
        auto itr_col( itr->Second().GetData().Begin());
        for( size_t col_number( 0); col_number < n_cols; ++col_number, ++itr_col)
        {
          contact_entropies( row_number)( CategoryEnum( table.GetHeader()( col_number))) = *itr_col;
        }
      }
      return contact_entropies;
    }
  } // namespace scorestat
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
// include header for this class
#include "scorestat/bcl_scorestat_protein_model_sse_triplet_chirality.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "assemble/bcl_assemble_voxel_grid_atom.h"
#include "biol/bcl_biol_aa_back_bone.h"
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
#include "score/bcl_score_aa_pair_hi_res_clash.h"
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
    const util::SiPtr< const util::ObjectInterface> ProteinModelSSETripletChirality::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new ProteinModelSSETripletChirality())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    ProteinModelSSETripletChirality *ProteinModelSSETripletChirality::Clone() const
    {
      return new ProteinModelSSETripletChirality( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelSSETripletChirality::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &ProteinModelSSETripletChirality::GetOutFilePostfix() const
    {
      static const std::string s_name( "sse_triplet_chirality"), s_con_name( "sse_triplet_chirality_contact"), s_rows_name( "sse_triplet_chirality_per_protein");
      return m_OutputPerProtein ? s_rows_name : ( m_ConsiderWhichSSEsInContact ? s_con_name : s_name);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ProteinModelSSETripletChirality::GetAlias() const
    {
      static const std::string s_name( "SSETripletChirality");
      return s_name;
    }

    //! @brief helper function to cache orientational information
    //! @param CACHE the cache to retrieve/store data in
    //! @param SSE_A_ID index of SSE_A
    //! @param SSE_B_ID index of SSE_B
    //! @param SSE_A, SSE_B the two sses of interest
    //! @return the orientation
    const assemble::SSEGeometryPacking::OrientationEnum &ProteinModelSSETripletChirality::GetCacheOrientation
    (
      storage::Vector< storage::Vector< assemble::SSEGeometryPacking::OrientationEnum> > &CACHE,
      const size_t &SSE_A_ID,
      const size_t &SSE_B_ID,
      const assemble::SSE &SSE_A,
      const assemble::SSE &SSE_B
    ) const
    {
      assemble::SSEGeometryPacking::OrientationEnum &ori( CACHE( SSE_A_ID)( SSE_B_ID));
      if( ori == assemble::SSEGeometryPacking::s_NumberOrientations)
      {
        ori = assemble::SSEGeometryPacking::OrientationEnum
              (
                assemble::SSEGeometryPacking::OrientationFromSSEs( SSE_A, SSE_B)
              );
      }
      return ori;
    }

    //! @brief hash a given packing. The given string is intended to be fast to hash, not necessarily easily readable
    std::string ProteinModelSSETripletChirality::GetPackingTripletHash
    (
      const biol::SSType &TYPE_A,
      const biol::SSType &TYPE_B,
      const biol::SSType &TYPE_C,
      const assemble::SSEGeometryPacking::OrientationEnum &PACKING_AB,
      const assemble::SSEGeometryPacking::OrientationEnum &PACKING_BC,
      const assemble::SSEGeometryPacking::OrientationEnum &PACKING_AC,
      const size_t &COUNTS_PACK_A,
      const size_t &COUNTS_PACK_B,
      const size_t &COUNTS_PACK_C,
      const bool   &ADJACENT,
      const bool   &RHS
    ) const
    {
      const size_t min_atoms_a
      (
        TYPE_A == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
      );
      const size_t min_atoms_b
      (
        TYPE_B == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
      );
      const size_t min_atoms_c
      (
        TYPE_C == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
      );
      return std::string
      (
        std::string( size_t( 1), TYPE_A.GetName()[ 0])
        + std::string( size_t( 1), TYPE_B.GetName()[ 0])
        + std::string( size_t( 1), TYPE_C.GetName()[ 0])
        + std::string( size_t( 1), PACKING_AB.GetString()[ 0])
        + std::string( size_t( 1), PACKING_BC.GetString()[ 0])
        + std::string( size_t( 1), PACKING_AC.GetString()[ 0])
        + (
            m_ConsiderWhichSSEsInContact
            ? util::Format()( COUNTS_PACK_A >= min_atoms_a * min_atoms_b)
              + util::Format()( COUNTS_PACK_B >= min_atoms_b * min_atoms_c)
              + util::Format()( COUNTS_PACK_C >= min_atoms_a * min_atoms_c)
            : std::string()
          )
        + std::string( ADJACENT ? "Adj" : "NAdj")
        + std::string( size_t( 1), RHS ? 'R' : 'L')
      );
    }

    //! @brief hash a given packing into a size_t for speed. The size_t will have the range 0 - 223 ( 7 * 8 * 4)
    size_t ProteinModelSSETripletChirality::GetPackingTripletNumber
    (
      const biol::SSType &TYPE_A,
      const biol::SSType &TYPE_B,
      const biol::SSType &TYPE_C,
      const assemble::SSEGeometryPacking::OrientationEnum &PACKING_AB,
      const assemble::SSEGeometryPacking::OrientationEnum &PACKING_BC,
      const assemble::SSEGeometryPacking::OrientationEnum &PACKING_AC,
      const size_t &COUNTS_PACK_A,
      const size_t &COUNTS_PACK_B,
      const size_t &COUNTS_PACK_C,
      const bool   &ADJACENT,
      const bool   &RHS
    ) const
    {
      const size_t min_atoms_a
      (
        TYPE_A == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
      );
      const size_t min_atoms_b
      (
        TYPE_B == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
      );
      const size_t min_atoms_c
      (
        TYPE_C == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
      );
      if( m_ConsiderWhichSSEsInContact)
      {
        bool ncontact_ab( COUNTS_PACK_A < min_atoms_a * min_atoms_b);
        bool ncontact_ac( COUNTS_PACK_C < min_atoms_a * min_atoms_c);
        bool ncontact_bc( COUNTS_PACK_B < min_atoms_b * min_atoms_c);
        if( ( ncontact_ab && ( ncontact_ac || ncontact_bc)) || ( ncontact_ac && ncontact_bc))
        {
          ncontact_ab = !ncontact_ab;
          ncontact_ac = !ncontact_ac;
          ncontact_bc = !ncontact_bc;
        }
        return
           ( TYPE_A == biol::GetSSTypes().STRAND ? 512 : 0)
         | ( TYPE_B == biol::GetSSTypes().STRAND ? 256 : 0)
         | ( TYPE_C == biol::GetSSTypes().STRAND ? 128 : 0)
         | ( PACKING_AB == assemble::SSEGeometryPacking::e_Parallel ? 64 : 0)
         | ( PACKING_BC == assemble::SSEGeometryPacking::e_Parallel ? 32 : 0)
         | ( PACKING_AC == assemble::SSEGeometryPacking::e_Parallel ? 16 : 0)
         | ( ADJACENT ? 8 : 0)
         | ( ncontact_ab ? 2 : ( ncontact_bc ? 4 : ( ncontact_ac ? 6 : 0)))
         | ( RHS ? 1 : 0);
      }
      return ( TYPE_A == biol::GetSSTypes().STRAND ? 128 : 0)
             | ( TYPE_B == biol::GetSSTypes().STRAND ? 64 : 0)
             | ( TYPE_C == biol::GetSSTypes().STRAND ? 32 : 0)
             | ( PACKING_AB == assemble::SSEGeometryPacking::e_Parallel ? 16 : 0)
             | ( PACKING_BC == assemble::SSEGeometryPacking::e_Parallel ? 8 : 0)
             | ( PACKING_AC == assemble::SSEGeometryPacking::e_Parallel ? 4 : 0)
             | ( ADJACENT ? 2 : 0)
             | ( RHS ? 1 : 0);
    }

    //! @brief convert a packing triplet number to a string
    std::string ProteinModelSSETripletChirality::GetPackingTripletString( const size_t &HASH_NUMBER) const
    {
      static storage::VectorND< 2, storage::Vector< std::string> > s_all_hashes;
      storage::Vector< std::string> &hashes( s_all_hashes( m_ConsiderWhichSSEsInContact ? 1 : 0));
      if( hashes.IsEmpty())
      {
        hashes.Resize( m_ConsiderWhichSSEsInContact ? s_NumberHashesContactSplit : s_NumberHashesNoContactSplit);
        const storage::Vector< biol::SSType> sstypes
        (
          storage::Vector< biol::SSType>::Create( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND)
        );
        const storage::Vector< assemble::SSEGeometryPacking::OrientationEnum> oris
        (
          storage::Vector< assemble::SSEGeometryPacking::OrientationEnum>::Create
          (
            assemble::SSEGeometryPacking::e_AntiParallel,
            assemble::SSEGeometryPacking::e_Parallel
          )
        );
        storage::Vector< storage::VectorND< 3, size_t> > counts
        (
          storage::Vector< storage::VectorND< 3, size_t> >::Create
          (
            storage::VectorND< 3, size_t>( 1000, 1000, 1000),
            storage::VectorND< 3, size_t>( 0, 1000, 1000),
            storage::VectorND< 3, size_t>( 1000, 0, 1000),
            storage::VectorND< 3, size_t>( 1000, 1000, 0)
          )
        );
        if( !m_ConsiderWhichSSEsInContact)
        {
          counts.Resize( 1);
        }
        // iterate through each allowed combination of ss types, orientation, adjaceny, and counts
        for( auto itr_a( sstypes.Begin()), ss_end( sstypes.End()); itr_a != ss_end; ++itr_a)
        {
          for( auto itr_b( sstypes.Begin()); itr_b != ss_end; ++itr_b)
          {
            for( auto itr_c( sstypes.Begin()); itr_c != ss_end; ++itr_c)
            {
              for( auto itr_d( oris.Begin()), ori_end( oris.End()); itr_d != ori_end; ++itr_d)
              {
                for( auto itr_e( oris.Begin()); itr_e != ori_end; ++itr_e)
                {
                  for( auto itr_f( oris.Begin()); itr_f != ori_end; ++itr_f)
                  {
                    for( auto itr_g( counts.Begin()), itr_g_end( counts.End()); itr_g != itr_g_end; ++itr_g)
                    {
                      for( int adjacent( 0); adjacent < 2; ++adjacent)
                      {
                        for( int rhs( 0); rhs < 2; ++rhs)
                        {
                          // associate the hash number with the string
                          hashes
                          (
                            GetPackingTripletNumber
                            (
                              *itr_a, *itr_b, *itr_c,
                              *itr_d, *itr_e, *itr_f,
                              itr_g->First(), itr_g->Second(), itr_g->Third(),
                              bool( adjacent),
                              bool( rhs)
                            )
                          ) = GetPackingTripletHash
                              (
                                *itr_a, *itr_b, *itr_c,
                                *itr_d, *itr_e, *itr_f,
                                itr_g->First(), itr_g->Second(), itr_g->Third(),
                                bool( adjacent),
                                bool( rhs)
                              );
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      return hashes( HASH_NUMBER);
    }

    //! @brief convert a packing triplet number to a string
    linal::Vector< double> ProteinModelSSETripletChirality::ComputeExpectedCounts
    (
      const linal::Matrix< double> &ADJ_CONTACT_SEEN, // Rows - ss triplet-type, Cols - contact type
      const linal::Matrix< double> &DIST_CONTACT_SEEN, // Rows - ss triplet-type, Cols - contact type
      const linal::Matrix< double> &PARALLEL_PROBS, // Rows - ss pair-type, Cols - SSE distance (1,2,3+)
      const linal::Vector< size_t> &SSE_TRIPLET_COUNTS_ADJACENT, // size - 7, sse triplet type, lowest value indicates first sse
      const linal::Vector< size_t> &SSE_TRIPLET_COUNTS // size - 7, sse triplet type, lowest value indicates first sse, 0 - helix, 1 - strand,
    ) const
    {
      linal::Vector< double> output( GetNumberHashes(), double( 0.0));
      const storage::Vector< biol::SSType> sstypes
      (
        storage::Vector< biol::SSType>::Create( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND)
      );
      const storage::Vector< assemble::SSEGeometryPacking::OrientationEnum> oris
      (
        storage::Vector< assemble::SSEGeometryPacking::OrientationEnum>::Create
        (
          assemble::SSEGeometryPacking::e_AntiParallel,
          assemble::SSEGeometryPacking::e_Parallel
        )
      );
      // compute probabilities for each orientation triplet. This arises from the integral
      // sin(phi)*(1/2-phi/(2*pi))/2 from 0 to pi/2 (in the case of odd # of parallel elements, the more common case)
      // and integral( sin(phi)*phi/(2*pi)/2 from 0 to pi/2 for the even # of parallel elements
      // So the index is the number of parallel sses mod two
      const double orientation_prob[ 2] = { 0.25 / math::g_Pi, 0.25 * ( math::g_Pi - 1.0) / math::g_Pi};

      storage::Vector< storage::VectorND< 3, size_t> > counts
      (
        storage::Vector< storage::VectorND< 3, size_t> >::Create
        (
          storage::VectorND< 3, size_t>( 1000, 1000, 1000),
          storage::VectorND< 3, size_t>( 0, 1000, 1000),
          storage::VectorND< 3, size_t>( 1000, 0, 1000),
          storage::VectorND< 3, size_t>( 1000, 1000, 0)
        )
      );
      if( !m_ConsiderWhichSSEsInContact)
      {
        counts.Resize( 1);
      }
      linal::Vector< double> expected_contact_probs( 8);
      linal::Vector< double> expected_ori_probs( 8);
      linal::Vector< double> relevant_contact_probs( counts.GetSize());
      for( int adjacent( 0); adjacent < 2; ++adjacent)
      {
        const linal::Matrix< double> &seen_matrix( adjacent ? ADJ_CONTACT_SEEN : DIST_CONTACT_SEEN);
        const int adj_index( adjacent ? 0 : 2);
        const int para_index( adjacent ? 1 : 2);
        const linal::Vector< size_t> &counts_vector( adjacent ? SSE_TRIPLET_COUNTS_ADJACENT : SSE_TRIPLET_COUNTS);
        // iterate through each allowed combination of ss types, orientation, adjaceny, and counts
        size_t int_type_a( 0);
        for( auto itr_a( sstypes.Begin()), ss_end( sstypes.End()); itr_a != ss_end; ++itr_a, ++int_type_a)
        {
          size_t int_type_b( 0);
          for( auto itr_b( sstypes.Begin()); itr_b != ss_end; ++itr_b, ++int_type_b)
          {
            const size_t int_type_ab( int_type_a + int_type_b * 2);
            size_t int_type_c( 0);
            for( auto itr_c( sstypes.Begin()); itr_c != ss_end; ++itr_c, ++int_type_c)
            {
              const size_t int_type_ac( int_type_a + int_type_c * 2);
              const size_t int_type_bc( int_type_b + int_type_c * 2);
              const size_t int_type_abc( int_type_a + 2 * ( int_type_b + 2 * int_type_c));
              expected_contact_probs = seen_matrix.GetRow( int_type_abc);
              expected_contact_probs.SetToSum( 1.0);
              for( size_t c( 0); c < size_t( 8); ++c)
              {;
                expected_ori_probs( c) =
                  ( c & 4 ? PARALLEL_PROBS( int_type_ab, adj_index) : 1.0 - PARALLEL_PROBS( int_type_ab, adj_index))
                  * ( c & 2 ? PARALLEL_PROBS( int_type_bc, adj_index) : 1.0 - PARALLEL_PROBS( int_type_bc, adj_index))
                  * ( c & 1 ? PARALLEL_PROBS( int_type_ac, para_index) : 1.0 - PARALLEL_PROBS( int_type_ac, para_index));
                expected_ori_probs( c) *= orientation_prob[ ( ( ( c & 4) >> 2) + ( ( c & 2) >> 1) + ( c & 1)) & 1];
              }
              expected_ori_probs.SetToSum( 1.0);
              if( !m_ConsiderWhichSSEsInContact)
              {
                relevant_contact_probs( 0) =
                  expected_contact_probs( 3) + expected_contact_probs( 5)
                  + expected_contact_probs( 6) + expected_contact_probs( 7);
              }
              else
              {
                relevant_contact_probs( 0) = expected_contact_probs( 7);
                relevant_contact_probs( 1) = expected_contact_probs( 3);
                relevant_contact_probs( 2) = expected_contact_probs( 5);
                relevant_contact_probs( 3) = expected_contact_probs( 6);
              }

              const double count( counts_vector( int_type_abc));
              size_t ori_index( 0);
              for( auto itr_d( oris.Begin()), ori_end( oris.End()); itr_d != ori_end; ++itr_d)
              {
                for( auto itr_e( oris.Begin()); itr_e != ori_end; ++itr_e)
                {
                  for( auto itr_f( oris.Begin()); itr_f != ori_end; ++itr_f, ++ori_index)
                  {
                    for
                    (
                      size_t count_vec_i( 0), n_counts( relevant_contact_probs.GetSize());
                      count_vec_i < n_counts;
                      ++count_vec_i
                    )
                    {
                      const storage::VectorND< 3, size_t> &cvec( counts( count_vec_i));
                      double expected_val( relevant_contact_probs( count_vec_i) * count / 2.0 * expected_ori_probs( ori_index));
                      output
                      (
                        GetPackingTripletNumber
                        (
                          *itr_a, *itr_b, *itr_c,
                          *itr_d, *itr_e, *itr_f,
                          cvec.First(), cvec.Second(), cvec.Third(),
                          bool( adjacent),
                          false
                        )
                      ) = expected_val;
                      output
                      (
                        GetPackingTripletNumber
                        (
                          *itr_a, *itr_b, *itr_c,
                          *itr_d, *itr_e, *itr_f,
                          cvec.First(), cvec.Second(), cvec.Third(),
                          bool( adjacent),
                          true
                        )
                      ) = expected_val;
                    }
                  }
                }
              }
            }
          }
        }
      }
      return output;
    }

  //////////////
  // operator //
  //////////////

    //! @brief add statistics about the given protein and chains to the internal table and histograms
    //! @param ENSEMBLE the protein ensemble the statistics of which to be computed
    std::string ProteinModelSSETripletChirality::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // write statistics
      std::ostringstream stream;

      // min probability to score two residues as being an almost-certain contact
      const double min_p( 0.04);

      const size_t n_hashes( m_ConsiderWhichSSEsInContact ? s_NumberHashesContactSplit : s_NumberHashesNoContactSplit);
      // initializes maps for storing histograms of distances between all types of amino acid pairs
      linal::Vector< double> triplet_packing_counts( n_hashes, 0.0);
      storage::Vector< size_t> triplet_packing_sse_counts( n_hashes, size_t( 0));
      storage::Vector< std::string> triplet_strings( n_hashes);
      linal::Vector< double> triplet_packing_counts_nic( n_hashes, 0.0);

      if( m_OutputPerProtein)
      {
        stream << "storage::Table<double>\tIsNative\t";
        for( size_t hash_id( 0); hash_id < n_hashes; ++hash_id)
        {
          triplet_strings( hash_id) = GetPackingTripletString( hash_id);
          stream << triplet_strings( hash_id) << '\t';
        }
        stream << '\n';
      }
      linal::Matrix< double> adj_seen( 8, 8, double( 0));
      linal::Matrix< double> dist_seen( 8, 8, double( 0));

      storage::Vector< math::RunningAverage< double> > adj_fract_parallel( 4);
      storage::Vector< math::RunningAverage< double> > semiadj_fract_parallel( 4);
      storage::Vector< math::RunningAverage< double> > dist_fract_parallel( 4);
      linal::Vector< size_t> total_sses_seen( 8, size_t( 0));
      linal::Vector< size_t> total_sses_seen_adj( 8, size_t( 0));

      storage::Set< biol::AtomType> types_of_interest( biol::GetAtomTypes().GetFirstSidechainAtomTypes());
      // Create a matrix that will hold the x,y, and z coordinates for terminii of the first & second strands
      linal::Matrix3x3< double> xyz_coordinates( 0.0); // make a matrix of size 3 X 3
      for( auto itr_ensemble( ENSEMBLE.Begin()), itr_ensemble_end( ENSEMBLE.End()); itr_ensemble != itr_ensemble_end; ++itr_ensemble)
      {
        // get current protein model
        util::ShPtr< assemble::ProteinModel> sp_current_model( *itr_ensemble);
        const assemble::ProteinModel &protein_model( *sp_current_model);

        // get current pdb name before all the loops start
        const util::ShPtr< util::Wrapper< std::string> > &model_name_ptr(
          protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");

        // instantiate sequences CaCb from pdb for each chain one
        BCL_MessageStd( "computing chirality and packing for " + model_name);

        for
        (
          auto itr_chain( ( *itr_ensemble)->GetChains().Begin()), itr_chain_end( ( *itr_ensemble)->GetChains().End());
          itr_chain != itr_chain_end;
          ++itr_chain
        )
        {
          // need at least three sses in the chain
          if( ( *itr_chain)->GetNumberSSEs() < 2)
          {
            continue;
          }
          // collect all non-coil sse's
          const util::SiPtrVector< const assemble::SSE> structured_sses
          (
            ( *itr_chain)->GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
          );

          if( structured_sses.GetSize() < 2)
          {
            continue;
          }
          storage::Vector< size_t> atoms_in_sse;
          for( auto itr_sse( structured_sses.Begin()), itr_sse_end( structured_sses.End()); itr_sse != itr_sse_end; ++itr_sse)
          {
            auto new_atoms( ( *itr_sse)->GetAtoms( types_of_interest));
            atoms_in_sse.PushBack( new_atoms.GetSize());
          }
          assemble::VoxelGridAA interactions_detector( m_InteractionDistance);
          linal::Matrix< float> interactions_matrix
          (
            interactions_detector.GetSSEInteractionMatrix
            (
              structured_sses,
              ( *itr_chain)->GetAminoAcids(),
              2,
              m_InteractionDistance,
              false,
              min_p,
              false
            )
          );
          // cache of the packing objects computed so far
          storage::Vector< storage::Vector< assemble::SSEGeometryPacking::OrientationEnum> > orientations
          (
            structured_sses.GetSize(),
            storage::Vector< assemble::SSEGeometryPacking::OrientationEnum>
            (
              structured_sses.GetSize(),
              assemble::SSEGeometryPacking::OrientationEnum( assemble::SSEGeometryPacking::s_NumberOrientations)
            )
          );

          // handle SSE-pair interaction-weight and parallel/anti-parallel bias
          size_t a( 0);
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              itr_a( structured_sses.Begin()),
              itr_end( structured_sses.End());
            itr_a != itr_end;
            ++itr_a, ++a
          )
          {
            const biol::SSType &type_a( ( *itr_a)->GetType());
            const size_t a_type( type_a == biol::GetSSTypes().STRAND ? 1 : 0);
            size_t b( a + 1);
            const size_t min_atoms_a
            (
              type_a == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
            );
            for
            (
              util::SiPtrVector< const assemble::SSE>::const_iterator itr_b( itr_a + 1);
              itr_b != itr_end;
              ++itr_b, ++b
            )
            {
              const size_t b_type( ( *itr_b)->GetType() == biol::GetSSTypes().STRAND ? 2 : 0);
              const biol::SSType &type_b( ( *itr_b)->GetType());
              const size_t min_atoms_b
              (
                type_b == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
              );

              const size_t interactions_ab( interactions_matrix( a, b));
              // for every pair of SSEs
              const assemble::SSEGeometryPacking::OrientationEnum &packing_ab
              (
                GetCacheOrientation( orientations, a, b, **itr_a, **itr_b)
              );

              const double is_parallel( packing_ab == assemble::SSEGeometryPacking::e_Parallel ? 1.0 : 0.0);
              if( itr_b == itr_a + 1)
              {
                adj_fract_parallel( b_type | a_type) += is_parallel;
              }
              else
              {
                if( itr_b == itr_a + 2)
                {
                  semiadj_fract_parallel( b_type | a_type) += is_parallel;
                }
                dist_fract_parallel( b_type | a_type) += is_parallel;
              }

              util::SiPtrVector< const assemble::SSE>::const_iterator itr_c( itr_b + 1);
              if( itr_c == itr_end)
              {
                continue;
              }

              size_t c( b + 1);

              for( ; itr_c != itr_end; ++itr_c, ++c)
              {
                const biol::SSType &type_c( ( *itr_c)->GetType());
                const size_t c_type( type_c == biol::GetSSTypes().STRAND ? 4 : 0);
                const size_t min_atoms_c
                (
                  type_c == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
                );
                const size_t interactions_ac( interactions_matrix( a, c));
                const size_t interactions_bc( interactions_matrix( b, c));
                const size_t adj_type
                (
                  ( interactions_ac < min_atoms_a * min_atoms_c ? 0 : 1)
                  + ( interactions_ab < min_atoms_a * min_atoms_b ? 0 : 4)
                  + ( interactions_bc < min_atoms_b * min_atoms_c ? 0 : 2)
                );
                // adjacent
                if( c == a + 2)
                {
                  ++total_sses_seen_adj( a_type | b_type | c_type);
                  adj_seen( a_type | b_type | c_type, adj_type) += 1.0;
                }
                else
                {
                  ++total_sses_seen( a_type | b_type | c_type);
                  dist_seen( a_type | b_type | c_type, adj_type) += 1.0;
                }

                const assemble::SSEGeometryPacking::OrientationEnum &packing_ac
                (
                  GetCacheOrientation( orientations, a, c, **itr_a, **itr_c)
                );
                const assemble::SSEGeometryPacking::OrientationEnum &packing_bc
                (
                  GetCacheOrientation( orientations, b, c, **itr_b, **itr_c)
                );
                const linal::Vector3D root_position( ( *itr_b)->GetCenter());

                xyz_coordinates.GetRow( 0).CopyValues( ( *itr_a)->GetMainAxis().GetStartPoint() - root_position);
                xyz_coordinates.GetRow( 1).CopyValues( ( *itr_a)->GetMainAxis().GetEndPoint() - root_position);
                xyz_coordinates.GetRow( 2).CopyValues( ( *itr_c)->GetCenter() - root_position);
                // Calculate the determinant of a 3 X 3 matrix that has rows sorted in descending order of priority.
                // Opposite orders will have opposite signs. The sign is assigned to a value of 1 and returned as the
                // value of the stereocenter.
                const float determinant( xyz_coordinates.Determinant());

                size_t triplet_number
                (
                  GetPackingTripletNumber
                  (
                    type_a,
                    type_b,
                    type_c,
                    packing_ab,
                    packing_bc,
                    packing_ac,
                    interactions_ab,
                    interactions_bc,
                    interactions_ac,
                    itr_a + 1 == itr_b && itr_b + 1 == itr_c,
                    determinant >= 0
                  )
                );

                if
                (
                  ( interactions_ac < min_atoms_a * min_atoms_c ? 1 : 0)
                  + ( interactions_ab < min_atoms_a * min_atoms_b ? 1 : 0)
                  + ( interactions_bc < min_atoms_b * min_atoms_c ? 1 : 0)
                  > 1
                )
                {
                  ++triplet_packing_counts_nic( triplet_number);
                  continue;
                }

                BCL_MessageVrb
                (
                  GetPackingTripletString( triplet_number)
                  + ( *itr_a)->GetIdentification() + " " + ( *itr_b)->GetIdentification() + " " + ( *itr_c)->GetIdentification()
                  + " " + util::Format()( interactions_ab)
                  + " " + util::Format()( interactions_ac)
                  + " " + util::Format()( interactions_bc)
                );
                triplet_packing_counts( triplet_number)
                    += double( interactions_ac + interactions_ab + interactions_bc) / double( atoms_in_sse( a) + atoms_in_sse( b) + atoms_in_sse( c));
                triplet_packing_sse_counts( triplet_number) += 1;
              }
            }
          }
        }
        if( m_OutputPerProtein && triplet_packing_counts.Sum())
        {
          stream << model_name << "\t1\t";
          for( size_t hash_id( 0); hash_id < n_hashes; ++hash_id)
          {
            stream << triplet_packing_counts( hash_id) << '\t';
          }
          stream << '\n';
          stream << model_name << "Mirror\t-1\t";
          for( size_t hash_id( 0); hash_id < n_hashes; ++hash_id)
          {
            stream << triplet_packing_counts( hash_id ^ size_t( 1)) << '\t';
          }
          stream << '\n';
          triplet_packing_counts = 0.0;
          triplet_packing_sse_counts.SetAllElements( 0);
        }
      }

      if( !m_OutputPerProtein)
      {
        stream << this->GetLabel().ToString() << '\n';
        stream << "# Triplet-Hash\tWeightedSSEs\tSSEs\tExpectedSSEsInContact\n";
        auto itr_sse( triplet_packing_sse_counts.Begin());
        size_t hash_id( 0);
        linal::Matrix< double> parallel_probs( 4, 3);
        for( int i( 0); i < 4; ++i)
        {
          parallel_probs( i, 0) = adj_fract_parallel( i).GetAverage();
          parallel_probs( i, 1) = semiadj_fract_parallel( i).GetAverage();
          parallel_probs( i, 2) = dist_fract_parallel( i).GetAverage();
        }

        linal::Vector< double> expected_counts
        (
          ComputeExpectedCounts( adj_seen, dist_seen, parallel_probs, total_sses_seen_adj, total_sses_seen)
        );
        for
        (
          auto itr( triplet_packing_counts.Begin()),
               itr_end( triplet_packing_counts.End());
          itr != itr_end;
          ++itr, ++itr_sse, ++hash_id
        )
        {
          stream << GetPackingTripletString( hash_id) << ' '
                 << *itr << ' '
                 << *itr_sse << ' '
                 << expected_counts( hash_id)
                 << '\n';
        }
      }

      return stream.str();

    } // end of operator

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelSSETripletChirality::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Computes chirality-based entropies");
      parameters.AddInitializer
      (
        "distance",
        "max distance between atoms to be considered in contact",
        io::Serialization::GetAgent( &m_InteractionDistance),
        util::Format()( ProteinModelSSETripletChirality().m_InteractionDistance)
      );
      parameters.AddInitializer
      (
        "min sc atoms strand",
        "minimum side chains atoms betwen contacting SSEs",
        io::Serialization::GetAgent( &m_MinAtomsInContactStrand),
        util::Format()( ProteinModelSSETripletChirality().m_MinAtomsInContactStrand)
      );
      parameters.AddInitializer
      (
        "min sc atoms helix",
        "minimum side chains atoms betwen contacting SSEs",
        io::Serialization::GetAgent( &m_MinAtomsInContactHelix),
        util::Format()( ProteinModelSSETripletChirality().m_MinAtomsInContactHelix)
      );
      parameters.AddInitializer
      (
        "consider which sses contact",
        "whether to consider which SSEs are in contact. This is rather important for several SSE types, but only rarely"
        "affects the preferred handedness",
        io::Serialization::GetAgent( &m_ConsiderWhichSSEsInContact),
        util::Format()( ProteinModelSSETripletChirality().m_ConsiderWhichSSEsInContact)
      );
      parameters.AddInitializer
      (
        "output per protein",
        "If true, output will be raw counts per protein, suitable for use in MinimizeScoreWeightset. Two rows will be "
        "output per protein model. One contains the real chirality, the other the mirror image. IsNative score will be "
        "output as 1 for the native model, -1 for the mirror",
        io::Serialization::GetAgent( &m_OutputPerProtein),
        util::Format()( ProteinModelSSETripletChirality().m_OutputPerProtein)
      );
      return parameters;
    }
  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_radius_of_gyration.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_radius_of_gyration.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> RadiusOfGyration::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new RadiusOfGyration())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &RadiusOfGyration::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
          "Table",
          "Histogram",
          GetStaticClassName< RadiusOfGyration::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &RadiusOfGyration::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extenstions[] =
      {
          "radius_of_gyration.tbl",
          "radius_of_gyration.histograms",
          GetStaticClassName< RadiusOfGyration::OutputOption>()
      };
      return s_output_filename_extenstions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RadiusOfGyration::RadiusOfGyration() :
        m_OutputOption( e_Table),
        m_ChainIds( "")
    {
      // nothing to do
    }

    //! @brief virtual copy constructor
    RadiusOfGyration *RadiusOfGyration::Clone() const
    {
      return new RadiusOfGyration( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &RadiusOfGyration::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &RadiusOfGyration::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &RadiusOfGyration::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RadiusOfGyration::GetAlias() const
    {
      static std::string s_Name( "RadiusOfGyration");
      return s_Name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string RadiusOfGyration::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // create table for holding radius of gyration statistics
      storage::Table< double> radius_gyration_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            "subunits", "nr_coordinates", "nr_aa", "rgyr_sqr", "density"
          )
        )
      );

      // iterate over the protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get the current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get current pdb name
        const util::ShPtr< util::Wrapper< std::string> >
          &sp_model_name( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( sp_model_name.IsDefined() ? io::File::RemovePath( sp_model_name->GetData()) : "");

        BCL_MessageStd( "Computing Rgyr terms on " + model_name);
        // get all chains of current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // create a vector of coordinates for defined first side chain atoms of all chains
        util::SiPtrVector< const linal::Vector3D> defined_first_sc_atom_coordinates_all_chains;

        // number of chains used in calculation
        size_t number_of_chains( 0);

        // total number of residues, even outside chains
        size_t number_res( 0);

        // iterate over chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
            chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains, if m_ChainIds us empty, then analyze all chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
             continue;
          }

          // get coordinate for all side chain atoms
          const util::SiPtrVector< const linal::Vector3D>
            all_first_sc_atom_coordinates_current_chain( ( *chain_itr)->GetAtomCoordinates( biol::GetAtomTypes().GetFirstSidechainAtomTypes()));

          if( all_first_sc_atom_coordinates_current_chain.IsEmpty())
          {
            continue;
          }

          auto sses( ( *chain_itr)->GetSSEs());
          number_res += sses.LastElement()->GetFirstAA()->GetSeqID() - sses.FirstElement()->GetFirstAA()->GetSeqID() + 1;

          // create a vector of coordinates for defined first side chain atoms on current chain
          util::SiPtrVector< const linal::Vector3D> defined_first_sc_atom_coordinates_current_chain;

          // push back the defined ones into the vector
          for
          (
            util::SiPtrVector< const linal::Vector3D>::const_iterator
              coordinate_itr( all_first_sc_atom_coordinates_current_chain.Begin()),
              coordinate_itr_end( all_first_sc_atom_coordinates_current_chain.End());
            coordinate_itr != coordinate_itr_end;
            ++coordinate_itr
          )
          {
            if( ( *coordinate_itr)->IsDefined())
            {
              defined_first_sc_atom_coordinates_current_chain.PushBack( *coordinate_itr);
            }
          }

          // add the vector for the current chain to the vector for all chains
          defined_first_sc_atom_coordinates_all_chains.Append( defined_first_sc_atom_coordinates_current_chain);

          // increment number of chains used in calculation
          ++number_of_chains;
        }

        // create variable for holding radius of gyration of current model
        double radius_of_gyration( 0.0);

        // if membrane flag is given
        if( biol::Membrane::GetFlagMembrane()->GetFlag())
        {
          // calculate the collapsed radius of gyration
          radius_of_gyration = score::RadiusOfGyration::SquareRadiusOfGyrationCollapsed
              (
                defined_first_sc_atom_coordinates_all_chains,
                current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
              );
        }
        // else soluble proteins
        else
        {
          // calculate radius of gyration in the standard way
          radius_of_gyration = coord::SquareRadiusOfGyration( defined_first_sc_atom_coordinates_all_chains);
        }

        const double atom_vdw_radius( 1.5);

        const double est_volume
        (
          std::max
          (
            std::max
            (
              coord::EstimateVolume( defined_first_sc_atom_coordinates_all_chains, 3.4, 10.0, atom_vdw_radius),
              1.0
            ),
            4.0 / 3.0 * math::g_Pi / 0.74 * number_res * atom_vdw_radius * atom_vdw_radius * atom_vdw_radius
          )
        );

        // insert into table
        radius_gyration_table.InsertRow
        (
          model_name,
          storage::Vector< double>::Create
          (
            number_of_chains,
            defined_first_sc_atom_coordinates_all_chains.GetSize(),
            number_res,
            radius_of_gyration,
            double( number_res) / est_volume
          )
        );
      }

      // write statistics
      std::ostringstream stream;

      // if output option is table
      if( m_OutputOption == e_Table)
      {
        radius_gyration_table.WriteFormatted( stream);
      }
      else if( m_OutputOption == e_Histogram)
      {
        // create histograms that holds radius of gyration statistics
        const storage::Map< std::string, math::Histogram> radius_of_gyration_histograms
        (
          score::RadiusOfGyration::HistogramsFromTable( radius_gyration_table)
        );

        // iterate over histograms
        for
        (
          storage::Map< std::string, math::Histogram>::const_iterator
            histogram_itr( radius_of_gyration_histograms.Begin()), histogram_itr_end( radius_of_gyration_histograms.End());
          histogram_itr != histogram_itr_end;
          ++histogram_itr
        )
        {
          stream << histogram_itr->first + " histogram" << '\n';
          stream << histogram_itr->second << '\n';
        }
      }

      return stream.str();
    } // end of operator ()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RadiusOfGyration::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes radius of gyration statistics."
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "Table"
      );
      return parameters;
    }
  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_sheet_template.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_chain.h"
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_collector_topology_interface.h"
#include "assemble/bcl_assemble_collector_topology_sheet.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_fold_template.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SheetTemplate::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new SheetTemplate())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &SheetTemplate::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "List",
        GetStaticClassName< SheetTemplate::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &SheetTemplate::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "sheet_template.list",
        GetStaticClassName< SheetTemplate::OutputOption>()
      };
      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SheetTemplate::SheetTemplate() :
      m_OutputOption( e_List),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    SheetTemplate *SheetTemplate::Clone() const
    {
      return new SheetTemplate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &SheetTemplate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &SheetTemplate::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &SheetTemplate::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SheetTemplate::GetAlias() const
    {
      static std::string s_name( "SheetTemplate");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string SheetTemplate::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // vector to hold sheet templates
      storage::Vector< assemble::FoldTemplate> sheet_templates;

      // initialize collector
      const util::ShPtr< assemble::CollectorTopologyInterface> sp_collector
      (
        new assemble::CollectorTopologySheet()
      );

      // iterate over the protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get pdb id
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");
        std::string pdb_id( io::File::RemovePath( model_name));
        std::transform( pdb_id.begin(), pdb_id.end(), pdb_id.begin(), toupper);

        // get all chains in current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // iterate over all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator
            chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // create a chain model from current chain
          const assemble::ProteinModel &chain_model( util::ShPtr< assemble::Chain>( ( *chain_itr)->HardCopy()));

          // collect sheets
          const util::ShPtrVector< assemble::Domain> &sheets( assemble::CollectorSheet().Collect( chain_model));

          BCL_MessageDbg( "number sheets found: " + util::Format()( sheets.GetSize()));
          size_t sheet_ctr( 0);

          // iterate over sheets
          for
          (
            util::ShPtrVector< assemble::Domain>::const_iterator
              sheet_itr( sheets.Begin()), sheet_itr_end( sheets.End());
            sheet_itr != sheet_itr_end;
            ++sheet_itr
          )
          {
            ++sheet_ctr;

            BCL_MessageDbg
            (
              "Sheet #" + util::Format()( sheet_ctr) + "\n"
              + ( *sheet_itr)->GetTopology()->GetOrderedIdentification()
            );

            // if less than 2 sses, skip
            if( ( *sheet_itr)->GetNumberSSEs() < 2)
            {
              continue;
            }

            // get the elements vector
            const util::SiPtrVector< const assemble::SSEGeometryInterface> elements_vector
            (
              ( *sheet_itr)->GetTopology()->GetElements()
            );

            // initialize ShPtrVector for geometries
            util::ShPtrVector< assemble::SSEGeometryPhiPsi> geometry_vector;

            // iterate over sses in current elements_vector
            for
            (
              util::SiPtrVector< const assemble::SSEGeometryInterface>::const_iterator
                element_itr( elements_vector.Begin()), element_itr_end( elements_vector.End());
              element_itr != element_itr_end;
              ++element_itr
            )
            {
              // get sses in current sheet
              const util::SiPtrVector< const assemble::SSE> &sheet_sses( ( *sheet_itr)->GetSSEs());

              // iterate over all sses in the sheet
              for
              (
                util::SiPtrVector< const assemble::SSE>::const_iterator
                  sheet_sse_itr( sheet_sses.Begin()), sheet_sse_itr_end( sheet_sses.End());
                sheet_sse_itr != sheet_sse_itr_end;
                ++sheet_sse_itr
              )
              {
                // if geometries are equal
                if( **element_itr == **sheet_sse_itr)
                {
                  geometry_vector.PushBack
                  (
                    util::ShPtr< assemble::SSEGeometryInterface>( new assemble::SSEGeometryPhiPsi( **sheet_sse_itr))
                  );
                }
              }
            }

            // create the fold template
            const assemble::FoldTemplate sheet_template
            (
              geometry_vector,
              sp_collector,
              pdb_id,
              false
            );

            // insert the template if it is defined
            if( sheet_template.HasDefinedGeometries())
            {
              sheet_templates.PushBack( sheet_template);
            }
          }
        } // end of iteration over all chains
      } // end of iteration over protein ensemble

      // write statistics
      std::stringstream stream;
      if( m_OutputOption == e_List)
      {
        stream << sheet_templates.GetSize() << '\n';
        // iterate over sheet_templates
        for
        (
          storage::Vector< assemble::FoldTemplate>::const_iterator
            sheet_template_itr( sheet_templates.Begin()), sheet_template_itr_end( sheet_templates.End());
          sheet_template_itr != sheet_template_itr_end;
          ++sheet_template_itr
        )
        {
          sheet_template_itr->WriteCompact( stream);
        }
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SheetTemplate::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes sheet template statistics."
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "List"
      );
      return parameters;
    }
  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_side_chain_distance.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SideChainDistance::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new SideChainDistance())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &SideChainDistance::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Histogram",
        GetStaticClassName< SideChainDistance::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &SideChainDistance::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "side_chain_distance.histogram",
        GetStaticClassName< SideChainDistance::OutputOption>()
      };
      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SideChainDistance::SideChainDistance() :
      m_OutputOption( e_Histogram),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    SideChainDistance *SideChainDistance::Clone() const
    {
      return new SideChainDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &SideChainDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &SideChainDistance::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &SideChainDistance::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SideChainDistance::GetAlias() const
    {
      static std::string s_name( "SideChainDistance");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string SideChainDistance::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // statistics for distance of side chain center of mass to sse main axis
      storage::Map< biol::SSType, storage::Map< biol::AAType, math::Histogram> > side_chain_distance;

      // initialize side_chain_distance
      for
      (
        biol::SSTypes::const_iterator
          sse_type_itr( biol::GetSSTypes().Begin()), sse_type_itr_end( biol::GetSSTypes().COIL.GetIterator());
        sse_type_itr != sse_type_itr_end;
        ++sse_type_itr
      )
      {
        // iterate over aa types
        for
        (
          biol::AATypes::const_iterator
            aa_type_itr( biol::GetAATypes().Begin()), aa_type_itr_end( biol::GetAATypes().End());
          aa_type_itr != aa_type_itr_end;
          ++aa_type_itr
        )
        {
          side_chain_distance[ *sse_type_itr][ *aa_type_itr] = math::Histogram( 0.0, 0.1, 100);
        }
      }

      // iterate over protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get protein pdb filename
        const util::ShPtr< util::Wrapper< std::string> >
          &sp_model_filename( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_filename( sp_model_filename.IsDefined() ? io::File::RemovePath( sp_model_filename->GetData()) : "");

        // get all chains of current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // iterate over all chains of current protein model
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator
            chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // get helices and strands the model
          const util::SiPtrVector< const assemble::SSE> &all_sses
          (
            current_protein_model.GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
          );

          // the pdb file might be missing sse entries
          if( all_sses.IsEmpty())
          {
            BCL_MessageStd( "Warning: " + model_filename + " probably does not have sse entries in the pdb file, skipping");
            continue;
          }

          // iterate over all sses
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            // skip sses that are shorter than minimum fragment length
            if( ( *sse_itr)->GetSize() <= ( *sse_itr)->GetType()->GetFragmentLength())
            {
              continue;
            }

            // create an iterator to the central amino acid residue
            biol::AASequence::const_iterator center_aa_itr
            (
              ( *sse_itr)->Begin() + ( ( ( *sse_itr)->GetType()->GetFragmentLength()) - 1) / 2
            );

            // create an iterator to the end amino acid residue
            biol::AASequence::const_iterator end_aa_itr( ( *sse_itr)->GetData().End());

            // iterate over fragments
            for
            (
              util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator
                frag_itr( ( *sse_itr)->GetFragments().Begin()), frag_itr_end( ( *sse_itr)->GetFragments().End());
              frag_itr != frag_itr_end && center_aa_itr != end_aa_itr;
              ++frag_itr, ++center_aa_itr
            )
            {
              // get the center of mass of the center amino acid residue
              const linal::Vector3D com_center_aa( ( *center_aa_itr)->CalculateCenterOfMassOfSideChain());

              // get the line segment of current fragment
              const coord::LineSegment3D lineseg_from_fragment( ( *frag_itr)->GetMainAxis());

              // get distance from center aa to sse fragment
              const double distance
              (
                coord::CalculateDistancePointFromLineSegment( lineseg_from_fragment, com_center_aa).First()
              );

              // insert
              side_chain_distance[ ( *sse_itr)->GetType()][ ( *center_aa_itr)->GetType()].PushBack( distance);
            }
          }
        } // end of iteration over all chains
      } // end of iteration over protein ensemble

      // write statistics
      std::stringstream stream;

      if( m_OutputOption == e_Histogram)
      {
        // iterate over sse type
        for
        (
          biol::SSTypes::const_iterator
            sse_type_itr( biol::GetSSTypes().Begin()), sse_type_itr_end( biol::GetSSTypes().COIL.GetIterator());
          sse_type_itr != sse_type_itr_end;
          ++sse_type_itr
        )
        {
          // write the name of sse type
          stream << ( *sse_type_itr)->GetName() << '\n';

          // iterate over aa types
          for
          (
            biol::AATypes::const_iterator
              aa_type_itr( biol::GetAATypes().Begin()), aa_type_itr_end( biol::GetAATypes().End());
            aa_type_itr != aa_type_itr_end;
            ++aa_type_itr
          )
          {
            // write name of amino acid residue
            stream << ( *aa_type_itr)->GetName() << '\n';

            // write histogram
            stream << side_chain_distance[ *sse_type_itr][ *aa_type_itr];
          }
        }
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SideChainDistance::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes statistics for side chain distance from main SSE axis."
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "Histogram"
      );

      return parameters;
    }
  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_sse_count.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSECount::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new SSECount())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &SSECount::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Table",
        GetStaticClassName< SSECount::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &SSECount::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extentions[] =
      {
        "sse_count.tbl",
        GetStaticClassName< SSECount::OutputOption>()
      };
      return s_output_filename_extentions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSECount::SSECount() :
      m_OutputOption( e_Table),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    SSECount *SSECount::Clone() const
    {
      return new SSECount( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &SSECount::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &SSECount::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &SSECount::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SSECount::GetAlias() const
    {
      static std::string s_Name( "SSECount");
      return s_Name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string SSECount::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // create table for holding sse count statistics
      storage::Table< double> sse_count_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            "HELIX", "STRAND", "COIL"
          )
        )
      );

      // statistics of sse count
      storage::Map< biol::SSType, storage::Map< biol::EnvironmentType, size_t> > sse_region_count;

      // iterate over all protein models in the ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        const assemble::ProteinModel &current_protein_model( ( **protein_model_itr));

        // get protein pdb filename
        const util::ShPtr< util::Wrapper< std::string> >
          &sp_model_filename( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_filename( sp_model_filename.IsDefined() ? io::File::RemovePath( sp_model_filename->GetData()) : "");

        // get membrane for current protein model
        const util::ShPtr< biol::Membrane> &sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        if( !sp_membrane.IsDefined())
        {
          BCL_MessageDbg( "Skip: " + util::Format()( model_filename) + " does not have an associated membrane");
          continue;
        }

        // get all chains in current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // iterate over all chains of current protein model
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator
            chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains, if m_ChainIds is empty, analyze all chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            BCL_MessageDbg( "Skip undesired chain: " + util::Format()( ( *chain_itr)->GetChainID()));
            continue;
          }

          // get all sse in current chain
          const util::SiPtrVector< const assemble::SSE> &all_sses( ( *chain_itr)->GetSSEs());

          // iterate over all sses in current chain
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            // get current sse type
            const biol::SSType &current_sse_type( ( *sse_itr)->GetType());

            // skip protein models without sse entries in the pdb file
            if( ( *chain_itr)->GetNumberSSEs() == 1 && current_sse_type == biol::GetSSTypes().COIL)
            {
              BCL_MessageStd( util::Format()( model_filename) + " probably does not have sse entries in the pdb file, skipping");
              continue;
            }

            // skip undefined sses
            if( !current_sse_type.IsDefined())
            {
              continue;
            }

            // get environment type of current sse
            const biol::EnvironmentType current_environment_type
            (
              sp_membrane.IsDefined() ?
                  sp_membrane->DetermineEnvironmentType( ( *sse_itr)->GetCenter()) :
                    biol::GetEnvironmentTypes().e_Solution
            );

            ++sse_region_count[ ( *sse_itr)->GetType()][ current_environment_type];
          } // end of iteration over all sses in current chain
        } // end of iteration over all chains in current protein model
      } // end of iteration over protein ensemble

      // store sse statistics to a table
      if( m_OutputOption == e_Table)
      {
        // iterate over environment types
        for
        (
          biol::EnvironmentTypes::const_iterator
            env_itr( biol::GetEnvironmentTypes().Begin()), env_itr_end( biol::GetEnvironmentTypes().End());
          env_itr != env_itr_end;
          ++env_itr
        )
        {
          sse_count_table.InsertRow
          (
            ( *env_itr)->GetName(),
            storage::Vector< double>::Create
            (
              sse_region_count[ biol::GetSSTypes().HELIX][ *env_itr],
              sse_region_count[ biol::GetSSTypes().STRAND][ *env_itr],
              sse_region_count[ biol::GetSSTypes().COIL][ *env_itr]
            )
          );
        }

      }

      // write statistics
      std::stringstream stream;

      if( m_OutputOption == e_Table)
      {
         sse_count_table.WriteFormatted( stream);
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSECount::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes statistics of occurrence of secondary structure elements."
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "Table"
      );

      return parameters;
    }
  } // namespace scorestat
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

// include header for this class
#include "scorestat/bcl_scorestat_sse_membrane_alignment.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_chain.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_data.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_interface.h"
#include "biol/bcl_biol_environment_types.h"
#include "biol/bcl_biol_membrane.h"
#include "biol/bcl_biol_ss_types.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "math/bcl_math_cubic_spline.h"
#include "math/bcl_math_histogram.h"
#include "score/bcl_score_sse_membrane_alignment.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEMembraneAlignment::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new SSEMembraneAlignment())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &SSEMembraneAlignment::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Histogram",
        GetStaticClassName< SSEMembraneAlignment::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &SSEMembraneAlignment::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extensions[] =
      {
        "sse_membrane_alignment.histograms",
        GetStaticClassName< SSEMembraneAlignment::OutputOption>()
      };
      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEMembraneAlignment::SSEMembraneAlignment() :
      m_OutputOption( e_Histogram),
      m_NumberOfBins( 9),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    SSEMembraneAlignment *SSEMembraneAlignment::Clone() const
    {
      return new SSEMembraneAlignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &SSEMembraneAlignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &SSEMembraneAlignment::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns the number of bins for the histogram
    //! @return the number of bins for the histogram
    const size_t &SSEMembraneAlignment::GetNumberOfBins() const
    {
      return m_NumberOfBins;
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &SSEMembraneAlignment::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SSEMembraneAlignment::GetAlias() const
    {
      static std::string s_name( "SSEMembraneAlignment");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string SSEMembraneAlignment::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // statistics of sse membrane alignment
      storage::Map< biol::SSType, storage::Map< biol::EnvironmentType, math::Histogram> > sse_membrane_alignment;
      storage::Map< biol::EnvironmentType, storage::Vector< math::Histogram> > strand_axis_alignment;

      // iterate over SSTypes of interest
      for
      (
        biol::SSTypes::const_iterator
          ss_type_itr( biol::GetSSTypes().Begin()),
          ss_type_itr_end( biol::GetSSTypes().COIL.GetIterator());
        ss_type_itr != ss_type_itr_end;
        ++ss_type_itr
      )
      {
        // iterate over EnvironmentTypes
        for
        (
          biol::EnvironmentTypes::const_iterator
            env_type_itr( biol::GetEnvironmentTypes().Begin()), env_type_itr_end( biol::GetEnvironmentTypes().End());
          env_type_itr != env_type_itr_end;
          ++env_type_itr
        )
        {
          // initialize sse_membrane_alignment
          sse_membrane_alignment[ *ss_type_itr][ *env_type_itr] = math::Histogram
          (
            0,
            math::g_Pi / ( m_NumberOfBins * 2),
            m_NumberOfBins
          );

          // if the SSEType is Strand
          if( *ss_type_itr == biol::GetSSTypes().STRAND)
          {
            // initialize strand_axis_alignment
            strand_axis_alignment[ *env_type_itr] = storage::Vector< math::Histogram>
            (
              coord::GetAxes().GetEnumCount(),
              math::Histogram( 0, math::g_Pi / ( m_NumberOfBins * 2), m_NumberOfBins)
            );
          }
        } // end of iteration over EnvironmentTypes
      } // end of iteration over SSTypes

      // iterate over the protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get the current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get current pdb name before all loops start
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");

        // get the membrane associated with the current protein model
        const util::ShPtr< biol::Membrane> sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        // get all chains in current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // iterate over all chains in current protein model
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator
            chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            BCL_MessageStd( "Skip chain " + util::Format()( ( *chain_itr)->GetChainID()) + " in " + model_name);
            continue;
          }

          // get all sses in current chain
          const util::SiPtrVector< const assemble::SSE> &all_sses( ( *chain_itr)->GetSSEs());

          // iterate over all sses in current chain
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            // get all sse fragments of current sse
            const util::SiPtrVector< const assemble::SSEGeometryInterface> &all_fragments( ( *sse_itr)->GetSSEGeometries());

            // iterate over fragments of current sse
            for
            (
              util::SiPtrVector< const assemble::SSEGeometryInterface>::const_iterator
                fragment_itr( all_fragments.Begin()), fragment_itr_end( all_fragments.End());
              fragment_itr != fragment_itr_end;
              ++fragment_itr
            )
            {
              // get current sse
              const assemble::SSEGeometryInterface &current_sse( **fragment_itr);

              // determine the environment type in which the current sse can be found
              // get environment type of current sse
              const biol::EnvironmentType environment_type
              (
                sp_membrane.IsDefined() ?
                    sp_membrane->DetermineEnvironmentType( current_sse.GetCenter()) :
                      biol::GetEnvironmentTypes().e_Solution
              );

              if( !environment_type.IsDefined())
              {
                BCL_MessageCrt
                (
                  "environment for this sse is not defined " + current_sse.GetIdentification()
                );
                continue;
              }

              double align_weight( 1.0);
              score::SSEMembraneAlignment score_sse_membrane_alignment;

              // if the current sse is a strand
              if( current_sse.GetType() == biol::GetSSTypes().STRAND)
              {
                // change the align_weight for strand
                align_weight = score_sse_membrane_alignment.WeightXAxis( current_sse, *sp_membrane);

                // for strands the axes are important
                for
                (
                  coord::Axis::const_iterator
                    axis_itr( coord::GetAxes().Begin()), axis_itr_end( coord::GetAxes().End());
                  axis_itr != axis_itr_end;
                  ++axis_itr
                )
                {
                  strand_axis_alignment[ environment_type]( *axis_itr).PushBack
                  (
                    score_sse_membrane_alignment.AngleToMembranePlane( current_sse, *axis_itr, *sp_membrane)
                  );
                }
              }

              // insert the angle into the map with appropriate weight
              sse_membrane_alignment[ current_sse.GetType()][ environment_type].PushBack
              (
                score_sse_membrane_alignment.AngleToMembranePlane( current_sse, coord::GetAxes().e_Z, *sp_membrane),
                align_weight
              );
            } // end of itration over fragments of sse
          } // end of iteration over all sses in the current chain
        } // end of iteration over all chains in the current protein model
      } // end of iteration over the protein ensemble

      // write out the statistics
      std::stringstream stream;
      io::OFStream write;

      if( m_OutputOption == e_Histogram)
      {
        // write output files
        io::File::MustOpenOFStream( write, "sse_membrane_alignment.histograms");
        write << sse_membrane_alignment;
        io::File::CloseClearFStream( write);

        io::File::MustOpenOFStream( write, "strand_angle_membrane.histograms");
        write << strand_axis_alignment;
        io::File::CloseClearFStream( write);

        // print instructions
        stream << "Refer to the following two files: \n";
        stream << "sse_membrane_alignment.histograms" << '\n';
        stream << "strand_angle_membrane.histograms";
      }

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSEMembraneAlignment::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes secondary structure element (SSE) membrane alignment statistics."
      );

      parameters.AddInitializer
      (
        "number_bins",
        "number of bins for the histogram",
        io::Serialization::GetAgent( &m_NumberOfBins),
        "9"
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "Histogram"
      );

      return parameters;
    }
  } // namespace scorestat
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
// include header for this class
#include "scorestat/bcl_scorestat_sse_packing.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model_data.h"
#include "assemble/bcl_assemble_sse_geometry_packer_all_fragment_pairs.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram_2d.h"
#include "score/bcl_score_aa_pair_sidechain_interaction.h"
#include "score/bcl_score_loop.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEPacking::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new SSEPacking())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &SSEPacking::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
        "Table",
        "Histogram2D",
        GetStaticClassName< SSEPacking::OutputOption>()
      };

      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &SSEPacking::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_file_extentions[] =
      {
        "sse_packing.tbl",
        "sse_packing.histograms2D",
        GetStaticClassName< SSEPacking::OutputOption>()
      };

      return s_output_file_extentions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEPacking::SSEPacking() :
      m_OutputOption( e_Table),
      m_SSEDistanceBinSize( 1.0),
      m_SSEAngleNumberBins( 24),
      m_StrandDistanceBinSize( 0.25),
      m_StrandAngleNumberBins( 24),
      m_ChainIds( "")
    {
      // nothing to do
    }

    //! @brief virtual copy constructor
    SSEPacking *SSEPacking::Clone() const
    {
      return new SSEPacking( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &SSEPacking::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &SSEPacking::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns the sse distance bin size for the 2D histogram
    //! @return the sse distance bin size for the 2D histogram
    const double &SSEPacking::GetSSEDistanceBinSize() const
    {
      return m_SSEDistanceBinSize;
    }

    //! @brief returns the sse twist angle bin size for the 2D histogram
    //! @return the sse twist angle bin size for the 2D histogram
    const size_t &SSEPacking::GetSSEAngleNumberBins() const
    {
      return m_SSEAngleNumberBins;
    }

    //! @brief returns the strand distance bin size for the 2D histogram
    //! @return the strand distance bin size for the 2D histogram
    const double &SSEPacking::GetStrandDistanceBinSize() const
    {
      return m_StrandDistanceBinSize;
    }

    //! @brief returns the strand twist angle bin size for the 2D histogram
    //! @return the strand twist angle bin size for the 2D histogram
    const size_t &SSEPacking::GetStrandAngleNumberBins() const
    {
      return m_StrandAngleNumberBins;
    }

    //! @brief returns the fragment minimum interface length
    //! @return the fragment minimum interface length
    const double &SSEPacking::GetFragmentMinInterfaceLength() const
    {
      return m_FragmentMinInterfaceLength;
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &SSEPacking::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SSEPacking::GetAlias() const
    {
      static const std::string s_name( "SSEPacking");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string SSEPacking::operator ()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // initialize a table for storing sse packing angles and distances
      storage::Table< double> sse_packing_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            "sse_packing_type", "twist_angle", "shortest_distance"
          )
        )
      );

      // initialize a vector of 2D histograms for storing sse packing angles and distances
      storage::Vector< math::Histogram2D> sse_packing_angle_distance
      (
        contact::Types::s_NumberValidTypes,
        math::Histogram2D
        (
          storage::VectorND< 2, double>( -math::g_Pi, 0),
          storage::VectorND< 2, double>( 2 * math::g_Pi / m_SSEAngleNumberBins, m_SSEDistanceBinSize),
          storage::VectorND< 2, size_t>( m_SSEAngleNumberBins, size_t( 20 / m_SSEDistanceBinSize))
        )
      );

      // initialize a 2D histogram for storing strand-strand packing angle and distance
      math::Histogram2D strand_strand_packing_angle_distance
      (
        storage::VectorND< 2, double>( -math::g_Pi, 0),
        storage::VectorND< 2, double>( 2 * math::g_Pi / m_StrandAngleNumberBins, m_StrandDistanceBinSize),
        storage::VectorND< 2, size_t>( m_StrandAngleNumberBins, size_t( 20 / m_StrandDistanceBinSize))
      );

      // initialize a vector of 2D histograms for storing sse fragment packing angles
      storage::Vector< math::Histogram2D> sse_fragment_angle_distance
      (
        contact::Types::s_NumberValidTypes,
        math::Histogram2D
        (
          storage::VectorND< 2, double>( -math::g_Pi, 0),
          storage::VectorND< 2, double>( 2 * math::g_Pi / m_StrandAngleNumberBins, m_SSEDistanceBinSize),
          storage::VectorND< 2, size_t>( m_SSEAngleNumberBins, size_t( 20 / m_SSEDistanceBinSize))
        )
      );

      // initialize a 2D histogram for storing strand-strand fragment packing angle and distance
      math::Histogram2D strand_strand_fragment_angle_distance
      (
        storage::VectorND< 2, double>( -math::g_Pi, 0),
        storage::VectorND< 2, double>( 2 * math::g_Pi / m_StrandAngleNumberBins, m_StrandDistanceBinSize),
        storage::VectorND< 2, size_t>( m_SSEAngleNumberBins, size_t( 20 / m_StrandDistanceBinSize))
      );

      // iterate over protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        util::ShPtr< assemble::ProteinModel> sp_protein_model( *protein_model_itr);
        const assemble::ProteinModel &protein_model( *sp_protein_model);

        // get current pdb name before all loops start
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");

        // get all sses in the current model
        const util::SiPtrVector< const assemble::SSE> all_sses( protein_model.GetSSEs());

        // iterate over all sses
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator
            sse_a_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
          sse_a_itr != sse_itr_end;
          ++sse_a_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *sse_a_itr)->GetChainID()) == std::string::npos)
          {
            continue;
          }

          // skip coils
          if( ( *sse_a_itr)->GetType() == biol::GetSSTypes().COIL)
          {
            continue;
          }

          // get the second sse
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator sse_b_itr( sse_a_itr + 1);
             sse_b_itr != sse_itr_end;
            ++sse_b_itr
          )
          {
            // skip undesired chains
            if( !m_ChainIds.empty() && m_ChainIds.find( ( *sse_b_itr)->GetChainID()) == std::string::npos)
            {
              continue;
            }

            // skip coils
            if( ( *sse_b_itr)->GetType() == biol::GetSSTypes().COIL)
            {
              continue;
            }

            // skip sse pairs that are possibly from a single broken sse
            if( biol::CalculateSequenceDistance( **sse_a_itr, **sse_b_itr) <= 1)
            {
              continue;
            }

            // create an object of SSEGeometryPacking from sse_a and sse_b
            const assemble::SSEGeometryPacking sse_pack( **sse_a_itr, **sse_b_itr);

            // collect twist angle and distance and store them into corresponding histograms
            if( m_OutputOption == e_Histogram)
            {
              SSEPackingAngleDistance( sse_pack, sse_packing_angle_distance, strand_strand_packing_angle_distance);
            }
            else if( m_OutputOption == e_Table)
            {
              // collect twist angle and distance and store them into corresponding table
              sse_packing_table.InsertRow
              (
                model_name + ": " + ( *sse_a_itr)->GetChainID() + "_" + ( *sse_b_itr)->GetChainID(),
                storage::Vector< double>::Create
                (
                  sse_pack.GetContactType(),
                  sse_pack.GetTwistAngle(),
                  sse_pack.GetDistance()
                ),
                true
              );
            }

            // get the fragments for sse_a and sse_b
            storage::Vector< storage::List< assemble::SSEGeometryPacking> > fragment_packing_lists
            (
              assemble::SSEGeometryPackerAllFragmentPairs( m_FragmentMinInterfaceLength).operator ()
              (
                **sse_a_itr, **sse_b_itr
              )
            );

            // iterate over a list of fragment packing list
            for
            (
              storage::Vector< storage::List< assemble::SSEGeometryPacking> >::const_iterator
                fragment_packing_list_itr( fragment_packing_lists.Begin()), fragment_packing_list_itr_end( fragment_packing_lists.End());
              fragment_packing_list_itr != fragment_packing_list_itr_end;
              ++fragment_packing_list_itr
            )
            {
              // iterate over fragment packing list
              for
              (
                storage::List< assemble::SSEGeometryPacking>::const_iterator
                  fragment_packing_itr( fragment_packing_list_itr->Begin()), fragment_packing_itr_end( fragment_packing_list_itr->End());
                fragment_packing_itr != fragment_packing_itr_end;
                ++fragment_packing_itr
              )
              {
                if( m_OutputOption == e_Histogram)
                {
                  // collect twist angle and distance and store them into corresponding 2D histograms
                  SSEPackingAngleDistance( *fragment_packing_itr, sse_fragment_angle_distance, strand_strand_fragment_angle_distance);
                }
                else if( m_OutputOption == e_Table)
                {
                  // collect twist angle and distance and store them into corresponding table
                  sse_packing_table.InsertRow
                  (
                    model_name + ": " + ( *sse_a_itr)->GetChainID() + "_" + ( *sse_b_itr)->GetChainID(),
                    storage::Vector< double>::Create
                    (
                      ( *fragment_packing_itr).GetContactType(),
                      ( *fragment_packing_itr).GetTwistAngle(),
                      ( *fragment_packing_itr).GetDistance()
                    ),
                    true
                  );
                }
              }
            }
          } // end of iterating over the second sse
        } // end of iterating over all sses
      } // end of iterating over protein ensemble

      // write to files
      std::stringstream stream;
      io::OFStream write;

      if( m_OutputOption == e_Histogram)
      {
        // write sse_packing statistics
        io::File::MustOpenOFStream( write, "sse_angle_distance.histograms2D");
        write << assemble::SSEGeometryPacking::GetDefaultFragmentMinimalInterfaceLength() << '\n';
        write << contact::GetTypes().HELIX_HELIX << '\n';
        write << sse_packing_angle_distance( contact::GetTypes().HELIX_HELIX);
        write << contact::GetTypes().HELIX_SHEET << '\n';
        write << sse_packing_angle_distance( contact::GetTypes().HELIX_SHEET);
        write << contact::GetTypes().HELIX_STRAND << '\n';
        write << sse_packing_angle_distance( contact::GetTypes().HELIX_STRAND);
        write << contact::GetTypes().STRAND_STRAND << '\n';
        write << sse_packing_angle_distance( contact::GetTypes().STRAND_STRAND);
        write << contact::GetTypes().SHEET_SHEET << '\n';
        write << sse_packing_angle_distance( contact::GetTypes().SHEET_SHEET);
        io::File::CloseClearFStream( write);

        // write sse_angle_distance_strand_strand
        io::File::MustOpenOFStream( write, "strand_angle_distance.histograms2D");
        write << assemble::SSEGeometryPacking::GetDefaultFragmentMinimalInterfaceLength() << '\n';
        write << contact::GetTypes().STRAND_STRAND << '\n';
        write << strand_strand_packing_angle_distance << '\n';
        io::File::CloseClearFStream( write);

        // write sse_fragment statistics
        io::File::MustOpenOFStream( write, "sse_fragment_angle_distance.histograms2D");
        write << contact::GetTypes().HELIX_HELIX << "\n";
        write << sse_fragment_angle_distance( contact::GetTypes().HELIX_HELIX) << "\n";
        write << contact::GetTypes().HELIX_SHEET << "\n";
        write << sse_fragment_angle_distance( contact::GetTypes().HELIX_SHEET) << "\n";
        write << contact::GetTypes().HELIX_STRAND << "\n";
        write << sse_fragment_angle_distance( contact::GetTypes().HELIX_STRAND) << "\n";
        write << contact::GetTypes().STRAND_STRAND << "\n";
        write << sse_fragment_angle_distance( contact::GetTypes().STRAND_STRAND) << "\n";
        write << contact::GetTypes().SHEET_SHEET << "\n";
        write << sse_fragment_angle_distance( contact::GetTypes().SHEET_SHEET) << "\n";
        io::File::CloseClearFStream( write);

        // write strand_fragment packing statistics
        io::File::MustOpenOFStream( write, "strand_fragment_angle_distance.histograms2D");
        write << m_SSEDistanceBinSize << '\n';
        write << contact::GetTypes().STRAND_STRAND << '\n';
        write << strand_strand_fragment_angle_distance;
        io::File::CloseClearFStream( write);

        // print instructions
        stream << "Refer to the following four files: " << '\n';
        stream << "sse_angle_distance.histograms2D" << '\n';
        stream << "strand_angle_distance.histograms2D" << '\n';
        stream << "sse_fragment_angle_distance.histograms2D" << '\n';
        stream << "strand_fragment_angle_distance.histograms2D";
      }
      else if( m_OutputOption == e_Table)
      {
        // write everything into one file
        sse_packing_table.WriteFormatted( stream);
      }

      return stream.str();
    } // end of operator

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSEPacking::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes sse packing statistics."
      );

      parameters.AddInitializer
      (
        "sse_distance_bin_size",
        "the bin size of distance for the 2D histogram",
        io::Serialization::GetAgent( &m_SSEDistanceBinSize),
        "1.0"
      );

      parameters.AddInitializer
      (
        "sse_angle_number_bins",
        "number of bins for sse distance",
        io::Serialization::GetAgent( &m_SSEAngleNumberBins),
        "24"
      );

      parameters.AddInitializer
      (
        "strand_distance_bin_size",
        "the bin size of distance for the 2D histogram",
        io::Serialization::GetAgent( &m_StrandDistanceBinSize),
        "0.25"
      );

      parameters.AddInitializer
      (
        "strand_angle_number_bins",
        "the bin size of angle for the 2D histogram",
        io::Serialization::GetAgent( &m_StrandAngleNumberBins),
        "24"
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "Histogram2D"
      );

      return parameters;
    }

    //! @brief collects twist angle and shortest distance for sse packing types
    //! @param SSE_PACK the sse packing type
    //!        SSEPACKING_ANGLE_DISTANCE a vector of 2D histograms that store twist angle and shortest distance
    //!        of sse packing types
    //!        STRAND_STRAND_ANGLE_DISTANCE a 2D histogram that stores the twist angle and shortest distance of
    //!        strand_strand packing
    //! @return a vector of 2D histograms that store twist angle and shortest distance of sse packing types
    storage::Vector< math::Histogram2D> &SSEPacking::SSEPackingAngleDistance
    (
      const assemble::SSEGeometryPacking &SSE_PACK,
      storage::Vector< math::Histogram2D> &SSEPACKING_ANGLE_DISTANCE,
      math::Histogram2D &STRAND_STRAND_ANGLE_DISTANCE
    ) const
    {
      // make pair of twist angle and shortest distance
      storage::VectorND< 2, double> angle_distance( SSE_PACK.GetTwistAngle(), SSE_PACK.GetDistance());

      // get packing type
      contact::Type contact_type( SSE_PACK.GetContactType());

      // consider all sorts of packing types, stores twist angle and shortest distance into 2D histogram for each sse packing type
      // strand-only packing
      if
      (
        contact_type == contact::GetTypes().STRAND_STRAND
        || contact_type == contact::GetTypes().SHEET_SHEET
        || contact_type == contact::GetTypes().UNDEFINED_STRAND_STRAND
      )
      {
        // push angle_distance pair to the corresponding 2D histogram
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().STRAND_STRAND).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * SSE_PACK.GetStrandStrandPairingWeight()
        );
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().SHEET_SHEET).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * SSE_PACK.GetRelativePositionWeight()
        );
        STRAND_STRAND_ANGLE_DISTANCE.PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * SSE_PACK.GetStrandStrandPairingWeight()
        );
      }
      // helix-only packing
      else if( contact_type == contact::GetTypes().HELIX_HELIX)
      {
        // push angle_distance pair to the corresponding 2D histogram
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().HELIX_HELIX).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * SSE_PACK.GetRelativePositionWeight()
        );
      }
      // mixed packing
      else if
      (
        contact_type == contact::GetTypes().HELIX_STRAND
        || contact_type == contact::GetTypes().STRAND_HELIX
      )
      {
        // push angle_distance pair to the corresponding 2D histogram
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().HELIX_STRAND).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * ( 1 - SSE_PACK.GetRelativePositionWeight())
        );
      }
      else if
      (
        contact_type == contact::GetTypes().HELIX_SHEET
        || contact_type == contact::GetTypes().SHEET_HELIX
      )
      {
        // push angle_distance pair to the corresponding 2D histogram
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().HELIX_SHEET).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * SSE_PACK.GetRelativePositionWeight()
        );
      }
      // undefined packing
      else if
      (
        contact_type == contact::GetTypes().UNDEFINED_HELIX_STRAND
        || contact_type == contact::GetTypes().UNDEFINED_STRAND_HELIX
      )
      {
        // push angle_distance pair to the corresponding 2D histogram
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().HELIX_SHEET).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * SSE_PACK.GetRelativePositionWeight()
        );
        SSEPACKING_ANGLE_DISTANCE( contact::GetTypes().HELIX_STRAND).PushBack
        (
          angle_distance, SSE_PACK.GetInteractionWeight() * ( 1 - SSE_PACK.GetRelativePositionWeight())
        );
      }

      // return the vector of 2D histograms
      return SSEPACKING_ANGLE_DISTANCE;
    }
  } // namespace scorestat
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
// include header for this class
#include "scorestat/bcl_scorestat_sspred_agreement.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "math/bcl_math_piecewise_function.h"
#include "math/bcl_math_roc_curve.h"
#include "score/bcl_score_environment_predictions.h"
#include "score/bcl_score_sse_predictions.h"
#include "sspred/bcl_sspred_ci_phi_psi.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "sspred/bcl_sspred_methods.h"
#include "sspred/bcl_sspred_pdb.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSPredAgreement::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new SSPredAgreement())
    );

    //! @brief default constructor
    SSPredAgreement::SSPredAgreement() :
      m_Method( sspred::GetMethods().e_PSIPRED),
      m_ChainIds( "")
    {
      // nothing else to do
    }

    //! @brief virtual copy constructor
    SSPredAgreement *SSPredAgreement::Clone() const
    {
      return new SSPredAgreement( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &SSPredAgreement::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &SSPredAgreement::GetOutFilePostfix() const
    {
      static const std::string s_suffix( "sspred_agreement.histograms");
      return s_suffix;
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &SSPredAgreement::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SSPredAgreement::GetAlias() const
    {
      static std::string s_name( "SSPredAgreement");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string SSPredAgreement::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure that the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // declare and initialize containers for holding sspred agreement statistics
      // histograms for secondary structure confidences
      // stores for each sspred method, each native sstype, the probability histograms for all sstypes
      storage::Map< sspred::Method, storage::Map< biol::SSType, storage::List< storage::Pair< double, double> > > >
        sspred_native_roc;
      storage::Map< sspred::Method, storage::Map< biol::EnvironmentType, storage::List< storage::Pair< double, double> > > >
        env_native_roc;

      const storage::Vector< biol::EnvironmentType> &env_types( biol::GetEnvironmentTypes().GetReducedTypes());
      const storage::Vector< biol::SSType> &ss_types( biol::GetSSTypes().GetReducedTypes());

      // initialize these containers
      for
      (
        auto method_itr( m_Method.Begin()), method_itr_end( m_Method.End());
        method_itr != method_itr_end;
        ++method_itr
      )
      {
        auto &map_o_lists( sspred_native_roc[ *method_itr]);
        auto &map_o_env_lists( env_native_roc[ *method_itr]);

        // for all sstypes, put in seed values so that even if we don't see this state we can just assume that the PPV is =
        // the predictions
        for
        (
          auto ss_type_itr_a( ss_types.Begin()), ss_type_itr_a_end( ss_types.End());
          ss_type_itr_a != ss_type_itr_a_end;
          ++ss_type_itr_a
        )
        {
          map_o_lists[ *ss_type_itr_a].PushBack( storage::Pair< double, double>( 0.0, 0.0));
          map_o_lists[ *ss_type_itr_a].PushBack( storage::Pair< double, double>( 1.0, 1.0));
        }

        for
        (
          auto env_type_itr( env_types.Begin()), env_type_itr_end( env_types.End());
          env_type_itr != env_type_itr_end;
          ++env_type_itr
        )
        {
          map_o_env_lists[ *env_type_itr].PushBack( storage::Pair< double, double>( 0.0, 0.0));
          map_o_env_lists[ *env_type_itr].PushBack( storage::Pair< double, double>( 1.0, 1.0));
        }
      }

      // iterate over the protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get a modifiable copy of current protein model
        assemble::ProteinModel &current_protein_model( *( *protein_model_itr)->HardCopy());
        sspred::PDB::SetEnvironmentTypes( current_protein_model, true);
        sspred::CIPhiPsi().Calculate( current_protein_model, true);

        // get the membrane for current protein model
        const util::ShPtr< biol::Membrane> sp_membrane
        (
          current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        // get model name
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");
        const std::string final_pdb_path( io::File::SplitToPathAndFileName( model_name).First());
        const std::string pdb_id( io::File::RemoveLastExtension( io::File::RemovePath( model_name)));

        // get all chains in current protein model
        util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // iterate over all chains
        for
        (
          util::ShPtrVector< assemble::Chain>::iterator
            chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
            BCL_MessageStd( "Skip chain " + util::Format()( ( *chain_itr)->GetChainID()) + " in " + model_name);
            continue;
          }

          // create reference on current chain
          biol::AASequence &current_chain( *( *chain_itr)->GetSequence());

          // determine sspred methods for which prediction files are available
          const storage::Set< sspred::Method> available_methods
            ( sspred::MethodHandler::AvailablePredictionFiles( m_Method, ( *chain_itr)->GetChainID(), pdb_id, final_pdb_path));

          // print out message if there are prediction files missing
          if( m_Method.GetSize() != available_methods.GetSize())
          {
            BCL_MessageCrt
            (
              "not all requested SSPRED Methods have prediction files available for pdb and chain " +
              model_name + " " + ( *chain_itr)->GetChainID()
            );
          }

          // try to read in all available sspred files
          if( !sspred::MethodHandler::ReadPredictionsForAASequence( available_methods, current_chain, pdb_id, final_pdb_path))
          {
            BCL_MessageCrt
            (
              "can't read in prediction file for all requested SSPRED Methods for pdb and chain " +
              model_name + " " + ( *chain_itr)->GetChainID()
            );
          }
          else // if all files were successfully read in
          {
            // iterate over the amino acid residues
            for
            (
              biol::AASequence::const_iterator
                aa_itr( current_chain.Begin()), aa_itr_end( current_chain.End());
              aa_itr != aa_itr_end;
              ++aa_itr
            )
            {
              auto true_ss_env( ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_CIPhiPsi));
              if( !true_ss_env.IsDefined())
              {
                true_ss_env = ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_PDB);
              }
              if( !true_ss_env.IsDefined())
              {
                continue;
              }
              const biol::SSType native_ss_type( true_ss_env->GetOneStateSSPrediction());
              const biol::EnvironmentType native_env_type( true_ss_env->GetOneStateTMPrediction()->GetReducedType());
              for
              (
                auto method_itr( available_methods.Begin()), method_itr_end( available_methods.End());
                method_itr != method_itr_end;
                ++method_itr
              )
              {
                const util::SiPtr< const sspred::MethodInterface> &preds( ( *aa_itr)->GetSSPrediction( *method_itr));
                if( !preds.IsDefined())
                {
                  continue;
                }
                if( native_ss_type.IsDefined())
                {
                  auto sspreds( preds->GetThreeStatePrediction());
                  auto itr_sspreds( sspreds.Begin());
                  for( auto ss_itr( ss_types.Begin()), ss_itr_end( ss_types.End()); ss_itr != ss_itr_end; ++ss_itr, ++itr_sspreds)
                  {
                    sspred_native_roc[ *method_itr][ *ss_itr].PushBack( storage::Pair< double, double>( *itr_sspreds, *ss_itr == native_ss_type ? 1.0 : 0.0));
                  }
                }
                if( native_env_type.IsDefined())
                {
                  auto envs( preds->GetThreeStateTMPrediction());
                  auto itr_env( envs.Begin());
                  for
                  (
                    auto env_itr( env_types.Begin()), env_itr_end( env_types.End());
                    env_itr != env_itr_end;
                    ++env_itr, ++itr_env
                  )
                  {
                    env_native_roc[ *method_itr][ *env_itr].PushBack( storage::Pair< double, double>( *itr_env, *env_itr == native_env_type ? 1.0 : 0.0));
                  }
                }
              }
            }
          } // end of iteration over all chains
        } // end of iteration over current protein model
      } // end of iteration over protein ensemble

      storage::Map< sspred::Method, storage::Map< biol::SSType, math::PiecewiseFunction> >
        sspred_native_localppv;
      storage::Map< sspred::Method, storage::Map< biol::EnvironmentType, math::PiecewiseFunction> >
        env_native_localppv;

      // initialize these containers
      for
      (
        auto method_itr( m_Method.Begin()), method_itr_end( m_Method.End());
        method_itr != method_itr_end;
        ++method_itr
      )
      {
        auto &map_o_ss_lists( sspred_native_roc[ *method_itr]);
        auto &map_o_env_lists( env_native_roc[ *method_itr]);

        for( auto ss_itr( ss_types.Begin()), ss_itr_end( ss_types.End()); ss_itr != ss_itr_end; ++ss_itr)
        {
          sspred_native_localppv[ *method_itr][ *ss_itr] = math::ROCCurve( map_o_ss_lists[ *ss_itr], 0.5, true).GetLocalPPVCurve();
        }
        for( auto env_itr( env_types.Begin()), env_itr_end( env_types.End()); env_itr != env_itr_end; ++env_itr)
        {
          env_native_localppv[ *method_itr][ *env_itr] = math::ROCCurve( map_o_env_lists[ *env_itr], 0.5, true).GetLocalPPVCurve();
        }
      }

      // write statistics
      std::stringstream stream;

      // output file stream
      io::OFStream write;

      // write mean and sd over n residues
      io::File::MustOpenOFStream( write, score::SSEPredictions::GetDefaultHistogramFilename());
      write << sspred_native_localppv;
      io::File::CloseClearFStream( write);
      io::File::MustOpenOFStream( write, score::EnvironmentPredictions::GetDefaultHistogramFilename());
      write << env_native_localppv;
      io::File::CloseClearFStream( write);

      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSPredAgreement::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes statistics for secondary structure prediction agreement."
      );

      parameters.AddInitializer
      (
        "sspred_methods",
        "secondary structure prediction methods to be analyzed",
        io::Serialization::GetAgent( &m_Method),
        "(PSIPRED)"
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      return parameters;
    }

  } // namespace scorestat
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
// include header for this class
#include "scorestat/bcl_scorestat_strand_alignment.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model_data.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_packer_all_fragment_pairs.h"
#include "biol/bcl_biol_aa.h"
#include "biol/bcl_biol_aa_back_bone_completer.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram_2d.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "score/bcl_score_aa_pair_sidechain_interaction.h"
#include "score/bcl_score_loop.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> StrandAlignment::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new StrandAlignment())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &StrandAlignment::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_Names[] =
      {
        "Table",
        GetStaticClassName< StrandAlignment::OutputOption>()
      };

      return s_Names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &StrandAlignment::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_output_filename_extensions[] =
      {
        "strand_alignment.tbl",
        GetStaticClassName< StrandAlignment::OutputOption>()
      };

      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    StrandAlignment::StrandAlignment() :
        m_OutputOption( e_Table),
        m_VisualizationOutputPath( "./"),
        m_ChainIds( "")
    {
      // nothing to do
    }

    //! @brief virtual copy constructor
    StrandAlignment *StrandAlignment::Clone() const
    {
      return new StrandAlignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &StrandAlignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &StrandAlignment::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief gets the cutoff for distance between Carbonyl-Oxygen and Amide_Nitrogen
    //! @return the cutoff for distance between Carbonyl-Oxygen and Amide_Nitrogen
    const double &StrandAlignment::GetHydrogenBondOHCutoff() const
    {
      return m_HydrogenBondOHCutoff;
    }

    //! @brief gets the path where to write the pdbs to visualize measurements
    //! @return the path where to write the pdbs to visualize measurements
    const std::string &StrandAlignment::GetVisualizationOutputPath() const
    {
      return m_VisualizationOutputPath;
    }

    //! @brief gets chain ids
    //! @return chain ids
    const std::string &StrandAlignment::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief gets the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &StrandAlignment::GetAlias() const
    {
      const static std::string s_name( "StrandAlignment");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string StrandAlignment::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // table for output
      storage::Table< double> strand_pair_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            // topology=1,0 (anti)parallel; on=1,0 this line is a O-N distance and C-O-N angle;
            // oh=0,1 this line is a O-H distance and C-O-H angle, on+oh=1 always;
            "topology", "on", "oh", "seq_id_a", "seq_id_b", "angle_rad", "angle_degree", "cos_angle", "distance"
          )
        )
      );

      // iterate over protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current pdb name before all loops start
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( ( *protein_model_itr)->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");

        // make a copy of current protein model
        assemble::ProteinModel protein_model( *( *protein_model_itr)->Clone());

        // create complete model
        biol::AABackBoneCompleter backbone_completer( true, false, false);
        protein_model = *backbone_completer.CompleteProteinModel( protein_model);

        // get all sheets in the current model
        const assemble::CollectorSheet sheet_collector;
        const util::ShPtrVector< assemble::Domain> sheets( sheet_collector.Collect( protein_model));

        size_t sheet_number( 0);

        // iterate over all sheets
        for
        (
          util::ShPtrVector< assemble::Domain>::const_iterator
            sheet_itr( sheets.Begin()), sheet_itr_end( sheets.End());
          sheet_itr != sheet_itr_end;
          ++sheet_itr, ++sheet_number
        )
        {
          // get the current sheet
          const assemble::Domain &current_sheet( **sheet_itr);

          // use the graph representation of the sheet to find interacting strands
          const assemble::Topology::GraphType &current_sheet_graph( current_sheet.GetTopology()->GetGraph());

          // get a SiPtrVector of all sses in current_sheet
          const util::SiPtrVector< const assemble::SSE> current_sheet_sses( current_sheet.GetSSEs());

          // iterate over all sses in current_sheet
          size_t sse_a_number( 0);
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_itr_a( current_sheet_sses.Begin()), sse_itr_end( current_sheet_sses.End());
            sse_itr_a != sse_itr_end;
            ++sse_itr_a, ++sse_a_number
          )
          {
            size_t sse_b_number( sse_a_number + 1);

            // iterate over the rest sses
            for
            (
              util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr_b( sse_itr_a + 1);
                sse_itr_b != sse_itr_end;
              ++sse_itr_b, ++sse_b_number
            )
            {
              // skip strand pairs that are not interacting
              if
              (
                !current_sheet_graph.FindEdge
                (
                  *current_sheet_graph.FindVertex( *sse_itr_a), *current_sheet_graph.FindVertex( *sse_itr_b)
                ).IsDefined()
              )
              {
                continue;
              }

              // get the sse pair and packing information
              const assemble::SSE &strand_a( **sse_itr_a), &strand_b( **sse_itr_b);
              const assemble::SSEGeometryPacking strand_packing( strand_a, strand_b);
              const assemble::SSEGeometryPacking::OrientationEnum packing_orientation( strand_packing.GetOrientation());

              // skip undesired chains
              if( !m_ChainIds.empty() && m_ChainIds.find( strand_a.GetChainID()) == std::string::npos)
              {
                BCL_MessageVrb
                (
                  "Skip strand pair with chain IDs: " + util::Format()( strand_a.GetChainID())
                      + util::Format()( strand_b.GetChainID()) + ", the desired chains are: " + m_ChainIds
                );
                continue;
              }

              // skip unknown orientations
              if
              (
                packing_orientation != assemble::SSEGeometryPacking::e_Parallel
                && packing_orientation != assemble::SSEGeometryPacking::e_AntiParallel
              )
              {
                BCL_MessageVrb
                (
                  "Skip the unknown packing orientation: " + packing_orientation.GetString()
                );
                continue;
              }

              // separately (parallel and anti-parallel) collect atom pairs with potential hydrogen bonds
              std::multimap< int, int> parallel_pairs, antiparallel_pairs;

              // iterate over strand_a
              for
              (
                util::ShPtrVector< biol::AABase>::const_iterator
                  aa_itr_a( strand_a.Begin()), aa_itr_a_end( strand_a.End());
                aa_itr_a != aa_itr_a_end;
                ++aa_itr_a
              )
              {
                // iterate over strand_b
                for
                (
                  util::ShPtrVector< biol::AABase>::const_iterator
                    aa_itr_b( strand_b.Begin()), aa_itr_b_end( strand_b.End());
                  aa_itr_b != aa_itr_b_end;
                  ++aa_itr_b
                )
                {
                  // if strand_a and strand_b are anti-parallel
                  if( packing_orientation == assemble::SSEGeometryPacking::e_AntiParallel)
                  {
                    const int seq_id_a( ( **aa_itr_a).GetSeqID()), seq_id_b( ( **aa_itr_b).GetSeqID());
                    antiparallel_pairs.insert( std::make_pair( seq_id_a, seq_id_b));
                    antiparallel_pairs.insert( std::make_pair( seq_id_b, seq_id_a));
                  }
                  // else if strand_a and strand_b are parallel
                  else if( packing_orientation == assemble::SSEGeometryPacking::e_Parallel)
                  {
                    const int seq_id_a( ( **aa_itr_a).GetSeqID());
                    if( aa_itr_b != strand_b.Begin())
                    {
                      const int seq_id_b_prev( ( **( aa_itr_b - 1)).GetSeqID());
                      parallel_pairs.insert( std::make_pair( seq_id_b_prev, seq_id_a));
                    }
                    if( aa_itr_b != strand_b.End() - 1)
                    {
                      const int seq_id_b_next( ( **( aa_itr_b + 1)).GetSeqID());
                      parallel_pairs.insert( std::make_pair( seq_id_a, seq_id_b_next));
                    }
                  }
                } // end of strand_b
              } // end of stand_a

              // collect the returned data
              storage::List< storage::Vector< double> > data;

              // collect atoms to visualize
              storage::List< biol::Atom> atoms_to_visualize;

              typedef storage::List< storage::Pair< util::ShPtr< biol::AABase>, util::ShPtr< biol::AABase> > > HBList;

              // process anti-parallel strands
              HBList hbonds( SelectShortestHBonds( strand_a, strand_b, antiparallel_pairs));

              for
              (
                HBList::const_iterator hb_itr( hbonds.Begin()), hb_itr_end( hbonds.End());
                  hb_itr != hb_itr_end;
                ++hb_itr
              )
              {
                data.Append
                (
                  ComputeStrandData( hb_itr->First(), hb_itr->Second(), assemble::SSEGeometryPacking::e_AntiParallel)
                );
                atoms_to_visualize.Append( CollectStrandAtoms( hb_itr->First(), hb_itr->Second()));
              }

              // process parallel strands
              hbonds = SelectShortestHBonds( strand_a, strand_b, parallel_pairs);

              for
              (
                HBList::const_iterator hb_itr( hbonds.Begin()), hb_itr_end( hbonds.End());
                  hb_itr != hb_itr_end;
                ++hb_itr
              )
              {
                data.Append
                (
                  ComputeStrandData( hb_itr->First(), hb_itr->Second(), assemble::SSEGeometryPacking::e_Parallel)
                );
                atoms_to_visualize.Append( CollectStrandAtoms( hb_itr->First(), hb_itr->Second()));
              }

              // create the strand string
              const std::string strand_str
              (
                "sheet" + util::Format()( sheet_number) + "_" + packing_orientation.GetString()
                + "_strands" + util::Format()( sse_a_number) + "_" + util::Format()( sse_b_number)
              );

              // insert data into table
              if( m_OutputOption == e_Table)
              {
                for
                (
                  storage::List< storage::Vector< double> >::const_iterator
                    data_itr( data.Begin()), data_itr_end( data.End());
                  data_itr != data_itr_end;
                  ++data_itr
                )
                {
                  strand_pair_table.InsertRow( model_name + "_" + strand_str, *data_itr, true);
                }
              }

              // write pdb files with atoms to visualize
              if( !m_VisualizationOutputPath.empty())
              {
                // split model_name into path and file name
                storage::VectorND< 2, std::string> path_and_filename( io::File::SplitToPathAndFileName( model_name));

                // output_path + filename
                std::string output_path_and_filename
                (
                  m_VisualizationOutputPath + util::GetRuntimeEnvironment().GetPathSeperator()
                  + io::File::RemoveLastExtension( path_and_filename.Second()) + "_" + strand_str + "."
                  + io::File::GetLastExtension( path_and_filename.Second())
                );
                BCL_MessageStd( "Writing visualization file " + output_path_and_filename);

                // get all pdb lines from model
                pdb::Factory factory;
                util::ShPtrList< pdb::Line> visualization_file_lines( factory.WriteCompleteModelToPDBLines( protein_model));

                // add atom lines
                biol::AA amino_acid;
                std::size_t atom_number( visualization_file_lines.GetSize());
                for
                (
                  storage::List< biol::Atom>::const_iterator atom_itr( atoms_to_visualize.Begin()), atom_itr_end( atoms_to_visualize.End());
                    atom_itr != atom_itr_end;
                  ++atom_itr, ++atom_number
                )
                {
                  visualization_file_lines.Append( pdb::Factory::WriteAtomToLine( *atom_itr, amino_acid, 'Z', atom_number));
                }

                // create handler and add lines
                pdb::Handler pdb_handler;
                pdb_handler.AppendLines( visualization_file_lines);

                // write visualization pdb
                io::OFStream pdb_write_stream;
                io::File::MustOpenOFStream( pdb_write_stream, output_path_and_filename);
                pdb_handler.WriteLines( pdb_write_stream);
                io::File::CloseClearFStream( pdb_write_stream);

              } // end of writing pdb files for visualization
            } // end of iterating over the rest sses in current sheet
          } // end of iterating over all sses in current sheet
        } // end of iterating over sheets
      } // end of iterating over protein ensemble

      // write statistics
      std::ostringstream stream;
      if( m_OutputOption == e_Table)
      {
        strand_pair_table.WriteFormatted( stream);
      }
      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer StrandAlignment::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes strand alignment statistics."
      );

      parameters.AddInitializer
      (
        "hydrogen_bond_OH_cutoff",
        "maximum length of a hydrogen bond measured between O-H to be considered for pairing strands",
        io::Serialization::GetAgent( &m_HydrogenBondOHCutoff),
        "4.5"
      );

      parameters.AddInitializer
      (
        "write_strand_alignment_visualization_path",
        "path where to write pdb files that visualize the calculation for the strand alignment data",
        io::Serialization::GetAgent( &m_VisualizationOutputPath),
        "./"
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "Table"
      );

      return parameters;
    }

    //! @brief computes angle between the given three atoms; not using the linal functions b/c of Normalize() call
    //! @param ATOM_A the base atom for calculating vector A->B, A->C
    //! @param ATOM_B the first point
    //! @param ATOM_C the second point
    //! @return the angle in radians
    double StrandAlignment::ComputeAngle( const biol::Atom &ATOM_A, const biol::Atom &ATOM_B, const biol::Atom &ATOM_C) const
    {
      const linal::Vector3D vector_a_b( ( ATOM_B.GetCenter() - ATOM_A.GetCenter()).Normalize());
      const linal::Vector3D vector_a_c( ( ATOM_C.GetCenter() - ATOM_A.GetCenter()).Normalize());
      return linal::ProjAngle( vector_a_b, vector_a_c);
    }

    //! @brief computes the H-bond distance between C-O and H-N, returns undefined if H coordinates are not defined
    //! @brief this function use an oriented H-bond definition: C-O ---> H-N
    //! @param AA_A the first AA ShPtr, the O is taken from this
    //! @param AA_B the second AA ShPtr, the H is taken from this
    //! @return the distance; undefined if a shptr is undefined or if no coordinates for the H are present
    double StrandAlignment::ComputeHBondDistance( const util::ShPtr< biol::AABase> &AA_A, const util::ShPtr< biol::AABase> &AA_B) const
    {
      // if any of AA_A or AA_B is undefined, return nan
      if( !AA_A.IsDefined() || !AA_B.IsDefined())
      {
        return util::GetUndefinedDouble();
      }

      const biol::Atom &a_co( AA_A->GetAtom( biol::GetAtomTypes().O));
      const biol::Atom &b_nh( AA_B->GetAtom( biol::GetAtomTypes().H));

      // if the position of the hydrogen atom is not defined, return nan
      if( !b_nh.AllCoordinatesDefined())
      {
        return util::GetUndefinedDouble();
      }

      return linal::Distance( a_co.GetCenter(), b_nh.GetCenter());
    }

    //! @brief computes a list of values for a single H-bond between the two given AAs
    //! @brief H-bond are defined as oriented: C-O ---> H-N, so the C-O is from AA_A, the H-N is from AA_B
    //! @param AA_A the first AA ShPtr, the C-O is taken from this
    //! @param AA_B the second AA ShPtr, the H-N is taken from this
    //! @param ORIENTATION orientation of the sheet
    //! @return a list of vectors of H-bond specific values identical to the final table columns
    storage::List< storage::Vector< double> > StrandAlignment::ComputeStrandData
    (
      const util::ShPtr< biol::AABase> &AA_A,
      const util::ShPtr< biol::AABase> &AA_B,
      const assemble::SSEGeometryPacking::Orientation &ORIENTATION
    ) const
    {
      BCL_MessageDbg( "AA_A = " + util::Format()( AA_A.IsDefined()) + " AA_B = " + util::Format()( AA_B.IsDefined()));

      // a list of vectors to store strand data
      storage::List< storage::Vector< double> > strand_data;

      // skip undefined hydrogen bonds
      if( !AA_A.IsDefined() || !AA_B.IsDefined())
      {
        BCL_MessageStd( "Skip HBond: one or both of the AA ShPtr is not defined.");
        return strand_data;
      }

      // get interesting atoms C, O, H, N
      const biol::Atom aa_ac( AA_A->GetAtom( biol::GetAtomTypes().C));
      const biol::Atom aa_ao( AA_A->GetAtom( biol::GetAtomTypes().O));
      const biol::Atom aa_bh( AA_B->GetAtom( biol::GetAtomTypes().H));
      const biol::Atom aa_bn( AA_B->GetAtom( biol::GetAtomTypes().N));

      // skip undefined hydrogen atoms
      if( !aa_bh.AllCoordinatesDefined())
      {
        BCL_MessageStd( "Skip HBond: the hydrogen atom has no defined coordinates.");
        return strand_data;
      }

      // get "hydrogen bonding lengths"
      const double distance_ao_bh( linal::Distance( aa_ao.GetCenter(), aa_bh.GetCenter()));
      const double distance_ao_bn( linal::Distance( aa_ao.GetCenter(), aa_bn.GetCenter()));

      // only use O ---> H distance cutoff, more precise than O ---> N
      if( distance_ao_bh > m_HydrogenBondOHCutoff)
      {
        BCL_MessageVrb
        (
          "Skip HBond: O to H distance " + util::Format()( distance_ao_bh) + " > cutoff "
          + util::Format()( m_HydrogenBondOHCutoff)
        );
        return strand_data;
      }

      // get "hydrogen bond angles" and their cosine values
      const double angle_ac_ao_bn( ComputeAngle( aa_ao, aa_ac, aa_bn));
      const double cosine_ac_ao_bn( std::cos( angle_ac_ao_bn));
      BCL_MessageVrb
      (
        "HBond between AA SeqIDs " + util::Format()( AA_A->GetSeqID()) + " and " + util::Format()( AA_B->GetSeqID()) + "\n"
        "C--O--N angle = " + util::Format()( angle_ac_ao_bn) + " in radians or "
        + util::Format()( angle_ac_ao_bn / math::g_Pi * 180) + " in degrees \n"
        + " cosine = " + util::Format()( cosine_ac_ao_bn) + " O--N distance = " + util::Format()( distance_ao_bn) + "\n"
        + " orientation = " + util::Format()( ORIENTATION)
      );
      strand_data.PushBack
      (
        storage::Vector< double>::Create
        (
          ORIENTATION, 1, 0, AA_A->GetSeqID(), AA_B->GetSeqID(),
          angle_ac_ao_bn, angle_ac_ao_bn / math::g_Pi * 180, cosine_ac_ao_bn, distance_ao_bn
        )
      );

      const double angle_ac_ao_bh( ComputeAngle( aa_ao, aa_ac, aa_bh));
      const double cosine_ac_ao_bh( std::cos( angle_ac_ao_bh));

      BCL_MessageVrb
      (
        "HBond between AA SeqID " + util::Format()( AA_A->GetSeqID()) + " and " + util::Format()( AA_B->GetSeqID()) + "\n"
        "C--O--H angle = " + util::Format()( angle_ac_ao_bh) + " in radians or "
        + util::Format()( angle_ac_ao_bh / math::g_Pi * 180) + " in degrees \n"
        + " cosine = " + util::Format()( cosine_ac_ao_bh) + " O--H distance = " + util::Format()( distance_ao_bh) + "\n"
        + " orientation = " + util::Format()( ORIENTATION)
      );
      strand_data.PushBack
      (
        storage::Vector< double>::Create
        (
          ORIENTATION, 0, 1, AA_A->GetSeqID(), AA_B->GetSeqID(),
          angle_ac_ao_bh, angle_ac_ao_bh / math::g_Pi * 180, cosine_ac_ao_bh, distance_ao_bh
        )
      );

      return strand_data;
    } // end of ComputeStrandData

    //! @brief collect all atoms that are part of an H-bond; use an oriented H-bond definition: C-O ---> H-N
    //! @param AA_A the first AA ShPtr, the C-O is taken from this
    //! @param AA_B the second AA ShPtr, the H-N is taken from this
    //! @return the list of atoms of this H-bond
    storage::List< biol::Atom> StrandAlignment::CollectStrandAtoms
    (
      const util::ShPtr< biol::AABase> &AA_A, const util::ShPtr< biol::AABase> &AA_B
    ) const
    {
      BCL_MessageDbg( "AA_A = " + util::Format()( AA_A.IsDefined()) + " AA_B = " + util::Format()( AA_B.IsDefined()));

      // a list of atoms
      storage::List< biol::Atom> strand_atoms;

      // skip undefined amino acid residues
      if( !AA_A.IsDefined() || !AA_B.IsDefined())
      {
        return strand_atoms;
      }

      // get interesting atoms C, O, H, N
      const biol::Atom aa_ac( AA_A->GetAtom( biol::GetAtomTypes().C));
      const biol::Atom aa_ao( AA_A->GetAtom( biol::GetAtomTypes().O));
      const biol::Atom aa_bh( AA_B->GetAtom( biol::GetAtomTypes().H));
      const biol::Atom aa_bn( AA_B->GetAtom( biol::GetAtomTypes().N));

      // skip undefined hydrogen atoms
      if( !aa_bh.AllCoordinatesDefined())
      {
        BCL_MessageVrb( "Skip HBond: the hydrogen atom has no defined coordinates.");
        return strand_atoms;
      }

      // get "hydrogen bonding lengths"
      const double distance_ao_bh( biol::Distance( aa_ao, aa_bh));

      // only use O ---> H distance cutoff, more precise than O ---> N
      if( distance_ao_bh > m_HydrogenBondOHCutoff)
      {
        BCL_MessageVrb
        (
          "Skip HBond: O to H distance " + util::Format()( distance_ao_bh) + " > cutoff "
          + util::Format()( m_HydrogenBondOHCutoff)
        );
        return strand_atoms;
      }

      // push those atoms onto the list and return the list
      strand_atoms.PushBack( aa_ac);
      strand_atoms.PushBack( aa_ao);
      strand_atoms.PushBack( aa_bh);
      strand_atoms.PushBack( aa_bn);

      return strand_atoms;
    } // end of CollectStrandAtoms

    //! @brief selects the shortest possible H-bond (below a cutoff) for each AA from a given set of AA pairs
    //! @param STRAND_A first strand, one AA in each AA pair comes from this strand
    //! @param STRAND_B second strand, one AA in each AA pair comes from this strand
    //! @param POSSIBLE_HBONDS map of possible H-bonds (pairs of AA seqids)
    //! @return a list of pairs of AA ShPtr that are the shortest H-bonds below a cutoff
    storage::List< storage::Pair< util::ShPtr< biol::AABase>, util::ShPtr< biol::AABase> > > StrandAlignment::SelectShortestHBonds
    (
      const assemble::SSE &STRAND_A,
      const assemble::SSE &STRAND_B,
      const std::multimap< int, int> &POSSIBLE_HBONDS
    ) const
    {
      // a list of pairs of AA ShPtrs
      storage::List< storage::Pair< util::ShPtr< biol::AABase>, util::ShPtr< biol::AABase> > > shortest_hbonds;

      // iterate over SeqIDs for first amino acid residues
      for
      (
        std::multimap< int, int>::const_iterator
          hbonds_itr_a( POSSIBLE_HBONDS.begin()), hbonds_itr_end( POSSIBLE_HBONDS.end());
        hbonds_itr_a != hbonds_itr_end;
        hbonds_itr_a = POSSIBLE_HBONDS.upper_bound( hbonds_itr_a->first)
      )
      {
        // shortest distance for current key
        double current_key_shortest_distance( util::GetUndefinedDouble());

        // the ShPtr to AA with shortest distance to current key
        util::ShPtr< biol::AABase> current_key_shortest_distance_aa_sp;

        // get the iterator to the amino acid residue whose SeqID is current key, it could be from either STRAND_A or STRAND_B
        biol::AASequence::const_iterator first_aa_itr_a( STRAND_A.FindAABySeqID( hbonds_itr_a->first));
        biol::AASequence::const_iterator first_aa_itr_b( STRAND_B.FindAABySeqID( hbonds_itr_a->first));

        // this has to be &&, b/c the AA will only be in one of the strands, i.e. one of them always finds the end
        if( first_aa_itr_a == STRAND_A.End() && first_aa_itr_b == STRAND_B.End())
        {
          BCL_MessageStd( "Skipping pair, could not find first AA in SSEs")
          continue;
        }

        // get the ShPtr to the amino acid residue whose SeqID is current key
        util::ShPtr< biol::AABase> first_aa_sp( first_aa_itr_a != STRAND_A.End() ? *first_aa_itr_a : *first_aa_itr_b);

        if( !first_aa_sp.IsDefined())
        {
          BCL_MessageStd( "Skipping pair, first AA ShPtr is not defined")
          continue;
        }

        // iterate over SeqIDs for second amino acid residues
        for
        (
          std::multimap< int, int>::const_iterator
            hbonds_itr_b( hbonds_itr_a), hbonds_itr_b_end( POSSIBLE_HBONDS.upper_bound( hbonds_itr_a->first));
          hbonds_itr_b != hbonds_itr_b_end;
          ++hbonds_itr_b
        )
        {
          // get the iterator to the amino acid residue paired with the amino acid residue whose SeqID is current key
          // it could be from either STRAND_A or STRAND_B
          biol::AASequence::const_iterator second_aa_itr_a( STRAND_A.FindAABySeqID( hbonds_itr_b->second));
          biol::AASequence::const_iterator second_aa_itr_b( STRAND_B.FindAABySeqID( hbonds_itr_b->second));

          if( second_aa_itr_a == STRAND_A.End() && second_aa_itr_b == STRAND_B.End())
          {
            BCL_MessageVrb( "Skipping pair, could not find second AA in SSEs")
            continue;
          }

          // get the ShPtr to the amino acid residue paired with the amino acid residue whose SeqID is current key
          util::ShPtr< biol::AABase> second_aa_sp( second_aa_itr_a != STRAND_A.End() ? *second_aa_itr_a : *second_aa_itr_b);

          if( !second_aa_sp.IsDefined())
          {
            BCL_MessageVrb( "Skipping pair, second AA ShPtr is not defined")
            continue;
          }

          // get hydrogen bonding distance
          const double distance( ComputeHBondDistance( first_aa_sp, second_aa_sp));

          BCL_MessageVrb
          (
            "Distance (" + util::Format()( hbonds_itr_b->first)
            + ", " + util::Format()( hbonds_itr_b->second) + ") = " + util::Format()( distance)
          );

          // get the shortest distance
          if( distance < current_key_shortest_distance || !util::IsDefined( current_key_shortest_distance))
          {
            current_key_shortest_distance = distance;
            current_key_shortest_distance_aa_sp = second_aa_sp;
          }
        }

        if( !current_key_shortest_distance_aa_sp.IsDefined())
        {
          BCL_MessageVrb( "current_key_shortest_distance_aa_sp not defined");
          continue;
        }

        BCL_MessageVrb
        (
          "Shortest Distance (" + util::Format()( hbonds_itr_a->first)
          + ", " + util::Format()( current_key_shortest_distance_aa_sp->GetSeqID())
          + ") = " + util::Format()( current_key_shortest_distance)
        );

        shortest_hbonds.PushBack( storage::Pair< util::ShPtr< biol::AABase>, util::ShPtr< biol::AABase> >( first_aa_sp, current_key_shortest_distance_aa_sp));
      }

      return shortest_hbonds;
    } // end of SelectShortestHBonds

  } // namespace scorestat
} // namespace bcl

