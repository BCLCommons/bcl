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
#include "score/bcl_score_aa_pair_hi_res_clash.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_atom.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_histogram_2d.h"
#include "math/bcl_math_histogram_3d.h"
#include "score/bcl_score_energy_distribution.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_logger_interface.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    const util::SiPtr< const util::ObjectInterface> AAPairHiResClash::s_ClashInstance
    (
      util::Enumerated< AAPairDistanceInterface>::AddInstance( new AAPairHiResClash())
    );
    const util::SiPtr< const util::ObjectInterface> AAPairHiResClash::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new AAPairHiResClash())
    );

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &AAPairHiResClash::GetAlias() const
    {
      static const std::string s_name( "AAPairHiResClash");
      return s_name;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a specified histogram file
    //! @param SCHEME scheme to be used
    AAPairHiResClash::AAPairHiResClash() :
      m_HistogramFileName( "aapair_contact_probs.histograms3D"),
      m_Scheme( "aa_clash_hires"),
      m_Histograms(),
      m_AngularBinSize( 0.0),
      m_DistanceCutoff( 0.0),
      m_ConsiderLoops( true),
      m_InterfaceOnly( true)
    {
      // read the histogram file and store the energy functions
      ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AAPairHiResClash object that is copied from this one
    AAPairHiResClash *AAPairHiResClash::Clone() const
    {
      return new AAPairHiResClash( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAPairHiResClash::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate amino acid pairing potential for given amino acid pair
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @return amino acid pairing potential for given amino acid pair
    double AAPairHiResClash::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B
    ) const
    {
      // calculate the Cb distances
      const double cb_cb_distance( biol::FirstSidechainAtomDistance( AMINO_ACID_A, AMINO_ACID_B));

      // return the score
      return operator()( AMINO_ACID_A, AMINO_ACID_B, cb_cb_distance);
    }

    //! @brief Get the probability that two AAs have heavy atoms within 1A + VDW Radii of one another
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param DISTANCE distance between the amino acid pair
    //! @return amino acid pairing potential for given amino acid pair
    double AAPairHiResClash::GetContactProbability
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const double DISTANCE
    ) const
    {
      if( AMINO_ACID_A.GetType() > AMINO_ACID_B.GetType())
      {
        return GetContactProbability( AMINO_ACID_B, AMINO_ACID_A, DISTANCE);
      }

      // calculate the angles distance
      const linal::Vector3D &ca_a( AMINO_ACID_A.GetCA().GetCoordinates());
      const linal::Vector3D &cb_a( AMINO_ACID_A.GetFirstSidechainAtom().GetCoordinates());
      const linal::Vector3D &ca_b( AMINO_ACID_B.GetCA().GetCoordinates());
      const linal::Vector3D &cb_b( AMINO_ACID_B.GetFirstSidechainAtom().GetCoordinates());

      // check that distance is defined
      if
      (
        !util::IsDefined( DISTANCE)
        || !ca_a.IsDefined()
        || !ca_b.IsDefined()
        || !cb_a.IsDefined()
        || !cb_b.IsDefined()
        || !AMINO_ACID_A.GetType()->IsNaturalAminoAcid()
        || !AMINO_ACID_B.GetType()->IsNaturalAminoAcid()
      )
      {
        return 0.0;
      }

      const double angle_b1_a1_a2( linal::ProjAngleCosinus( ca_a, cb_a, cb_b));
      const double angle_b2_a2_a1( linal::ProjAngleCosinus( ca_b, cb_b, cb_a));

      double f
      (
        ( *m_Histograms)( AMINO_ACID_A.GetType())( AMINO_ACID_B.GetType())->Value( DISTANCE, angle_b2_a2_a1, angle_b1_a1_a2)
      );
      return std::max( f, 0.0);
    }

    //! @brief calculate amino acid pairing potential for given amino acid pair
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param DISTANCE distance between the amino acid pair
    //! @return amino acid pairing potential for given amino acid pair
    double AAPairHiResClash::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const double DISTANCE
    ) const
    {
      if( AMINO_ACID_A.GetType() > AMINO_ACID_B.GetType())
      {
        return operator()( AMINO_ACID_B, AMINO_ACID_A, DISTANCE);
      }

      // calculate the angles distance
      const linal::Vector3D &ca_a( AMINO_ACID_A.GetCA().GetCoordinates());
      const linal::Vector3D &cb_a( AMINO_ACID_A.GetFirstSidechainAtom().GetCoordinates());
      const linal::Vector3D &ca_b( AMINO_ACID_B.GetCA().GetCoordinates());
      const linal::Vector3D &cb_b( AMINO_ACID_B.GetFirstSidechainAtom().GetCoordinates());

      // check that distance is defined
      if
      (
        !util::IsDefined( DISTANCE)
        || !ca_a.IsDefined()
        || !ca_b.IsDefined()
        || !cb_a.IsDefined()
        || !cb_b.IsDefined()
        || !AMINO_ACID_A.GetType()->IsNaturalAminoAcid()
        || !AMINO_ACID_B.GetType()->IsNaturalAminoAcid()
      )
      {
        return 0.0;
      }

      const linal::Matrix< double> &clash_matrix( *( *m_ClashDistances)( AMINO_ACID_A.GetType())( AMINO_ACID_B.GetType()));
      const size_t nr_angular_bins( clash_matrix.GetNumberRows() - 1);
      const size_t angle_b1_a1_a2
      (
        std::min( size_t( ( linal::ProjAngleCosinus( ca_a, cb_a, cb_b) + 1.0) / m_AngularBinSize), nr_angular_bins)
      );
      const size_t angle_b2_a2_a1
      (
        std::min( size_t( ( linal::ProjAngleCosinus( ca_b, cb_b, cb_a) + 1.0) / m_AngularBinSize), nr_angular_bins)
      );

      const double clash_distance( clash_matrix( angle_b2_a2_a1, angle_b1_a1_a2));

      double f( DISTANCE < clash_distance ? clash_distance - DISTANCE : 0.0);
      const biol::AAType type_a( AMINO_ACID_A.GetType()), type_b( AMINO_ACID_B.GetType());
      for( auto itr_a( AMINO_ACID_A.GetAtoms().Begin()), itr_a_end( AMINO_ACID_A.GetAtoms().End()); itr_a != itr_a_end; ++itr_a)
      {
        if( !( *itr_a)->GetCoordinates().IsDefined() || !( *itr_a)->GetType()->IsBackBone())
        {
          continue;
        }
        const chemistry::ElementType h_bond_partner
        (
          ( *itr_a)->GetType()->GetElementType() == chemistry::GetElementTypes().e_Hydrogen
          ? chemistry::GetElementTypes().e_Oxygen
          : ( *itr_a)->GetType()->GetElementType() == chemistry::GetElementTypes().e_Oxygen
            ? chemistry::GetElementTypes().e_Hydrogen
            : chemistry::ElementType()
        );
        const double vdw_a( type_a->GetVdwRadiusToOtherAA( ( *itr_a)->GetType()));
        for( auto itr_b( AMINO_ACID_B.GetAtoms().Begin()), itr_b_end( AMINO_ACID_B.GetAtoms().End()); itr_b != itr_b_end; ++itr_b)
        {
          if( !( *itr_b)->GetCoordinates().IsDefined() || !( *itr_b)->GetType()->IsBackBone())
          {
            continue;
          }
          double dist( linal::Distance( ( *itr_b)->GetCoordinates(), ( *itr_a)->GetCoordinates()));
          // Hydrogen bonds can (rarely) have shorter distance than allowed by VdW radii. Handle them specially; just
          // ensure that the distance is less than a covalent bond between the atoms
          if( ( *itr_b)->GetType()->GetElementType() == h_bond_partner)
          {
            if( dist < 1.2)
            {
              f = std::max( f, 1.2 - dist);
            }
          }
          else
          {
            f = std::max( f, std::max( type_b->GetVdwRadiusToOtherAA( ( *itr_b)->GetType()) + vdw_a - dist, 0.0));
          }
        }
      }
      if( f > 0.0)
      {
        BCL_MessageVrb
        (
          "Clash between " + AMINO_ACID_A.GetIdentification()
          + " and " + AMINO_ACID_B.GetIdentification() + " " + util::Format()( f)
        );
      }
      return std::min( f, 1.0);
    }

    //! @brief calculate clash score for a protein model
    //! @param MODEL the protein model of interest
    //! @return amino acid pairing potential for given protein
    double AAPairHiResClash::operator()( const assemble::ProteinModel &MODEL) const
    {
      assemble::VoxelGridAA s_voxel_grid( m_DistanceCutoff);

      auto inter_sse_contacts
      (
        s_voxel_grid.GetSSEConnections
        (
          MODEL.GetSSEs(),
          MODEL.GetAminoAcids(),
          2,
          m_DistanceCutoff,
          m_ConsiderLoops
        )
      );
      double score( 0.0);
      for( auto itr( inter_sse_contacts.Begin()), itr_end( inter_sse_contacts.End()); itr != itr_end; ++itr)
      {
        const auto &triplet( *itr);
        score += this->operator()( *triplet.First(), *triplet.Second(), triplet.Third());
      }
      return score;
    }

    //! @brief calculate clash score for SSEs
    //! @param SSE_A, SSE_B the SSEs to check for clashes
    //! @return clash score
    double AAPairHiResClash::operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
    {
      if( m_InterfaceOnly && SSE_A.GetChainID() == SSE_B.GetChainID())
      {
        return 0.0;
      }
      assemble::VoxelGridAA voxel_grid_a( m_DistanceCutoff);
      voxel_grid_a.SetObjects( util::SiPtrVector< const biol::AABase>( SSE_A.Begin(), SSE_A.End()));
      double score( 0.0);
      for( auto itr_b( SSE_B.Begin()), itr_b_end( SSE_B.End()); itr_b != itr_b_end; ++itr_b)
      {
        if( !( *itr_b)->GetType()->IsNaturalAminoAcid())
        {
          continue;
        }
        storage::Vector< storage::Pair< util::SiPtr< const biol::AABase>, double> >
          res( voxel_grid_a.GetNeighbors( **itr_b, m_DistanceCutoff));
        for( auto itr( res.Begin()), itr_end( res.End()); itr != itr_end; ++itr)
        {
          const auto &pr( *itr);
          score += this->operator()( **itr_b, *pr.First(), pr.Second());
        }
      }
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    AAPairHiResClash::WriteDetailedSchemeAndValues
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      std::ostream &OSTREAM
    ) const
    {
      // calculate the CB Distance
      const double cb_distance( biol::FirstSidechainAtomDistance( AMINO_ACID_A, AMINO_ACID_B));

      // if distance is not defined return 0
      if( !util::IsDefined( cb_distance))
      {
        return OSTREAM;
      }
      // calculate the angles distance
      const linal::Vector3D &ca_a( AMINO_ACID_A.GetCA().GetCoordinates());
      const linal::Vector3D &cb_a( AMINO_ACID_A.GetFirstSidechainAtom().GetCoordinates());
      const linal::Vector3D &ca_b( AMINO_ACID_B.GetCA().GetCoordinates());
      const linal::Vector3D &cb_b( AMINO_ACID_B.GetFirstSidechainAtom().GetCoordinates());

      // check that distance is defined
      if
      (
        !ca_a.IsDefined()
        || !ca_b.IsDefined()
        || !cb_a.IsDefined()
        || !cb_b.IsDefined()
        || !AMINO_ACID_A.GetType()->IsNaturalAminoAcid()
        || !AMINO_ACID_B.GetType()->IsNaturalAminoAcid()
      )
      {
        return OSTREAM;
      }

      const double angle_b1_a1_a2( linal::ProjAngleCosinus( ca_a, cb_a, cb_b));
      const double angle_b2_a2_a1( linal::ProjAngleCosinus( ca_b, cb_b, cb_a));

      // write Scheme
      OSTREAM << AMINO_ACID_A.GetSeqID() << '\t'
              << AMINO_ACID_A.GetType()->GetThreeLetterCode() << '\t'
              << AMINO_ACID_B.GetSeqID() << '\t'
              << AMINO_ACID_B.GetType()->GetThreeLetterCode() << '\t'
              << cb_distance << '\t'
              << angle_b2_a2_a1 << '\t'
              << angle_b1_a1_a2 << '\t'
              << operator()( AMINO_ACID_A, AMINO_ACID_B, cb_distance) << '\n';

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param MODEL model of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    AAPairHiResClash::WriteDetailedSchemeAndValues
    (
      const assemble::ProteinModel &MODEL,
      std::ostream &OSTREAM
    ) const
    {
      assemble::VoxelGridAA s_voxel_grid( m_DistanceCutoff);
      auto inter_sse_contacts
      (
        s_voxel_grid.GetSSEConnections
        (
          MODEL.GetSSEs(),
          MODEL.GetAminoAcids(),
          2,
          m_DistanceCutoff,
          false
        )
      );
      for( auto itr( inter_sse_contacts.Begin()), itr_end( inter_sse_contacts.End()); itr != itr_end; ++itr)
      {
        const auto &triplet( *itr);
        this->WriteDetailedSchemeAndValues( *triplet.First(), *triplet.Second(), OSTREAM);
      }
      return OSTREAM;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAPairHiResClash::GetSerializer() const
    {
      io::Serializer serializer;

      serializer.SetClassDescription
      (
        "Scores amino acid pair distance and angles for clash "
      );
      serializer.AddInitializer
      (
        "consider loops",
        "whether consider loops when computing contact energy and clashes",
        io::Serialization::GetAgent( &m_ConsiderLoops),
        "True"
      );
      serializer.AddInitializer
      (
        "interface only",
        "whether to score interface only",
        io::Serialization::GetAgent( &m_InterfaceOnly),
        "True"
      );

      return serializer;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    namespace
    {
      //! @brief get a cached clash distance map
      //! @param HISTOGRAMS the histogram3ds to use to derive the clash distances
      const storage::Vector< util::ShPtrVector< linal::Matrix< double> > > &GetClashDistances
      (
        const util::SiPtr< const storage::Vector< util::ShPtrVector< math::Histogram3D> > > &HISTOGRAM3DS
      )
      {
        static storage::Map
        <
          util::SiPtr< const storage::Vector< util::ShPtrVector< math::Histogram3D> > >,
          storage::Vector< util::ShPtrVector< linal::Matrix< double> > >
        > s_histograms;
        storage::Vector< util::ShPtrVector< linal::Matrix< double> > > &histogram_vec( s_histograms[ HISTOGRAM3DS]);
        if( !histogram_vec.IsEmpty())
        {
          return histogram_vec;
        }
        // store map of all aa pairs with their histogram
        histogram_vec.Reset();
        histogram_vec.Resize
        (
          biol::AATypes::s_NumberStandardAATypes,
          util::ShPtrVector< linal::Matrix< double> >( biol::AATypes::s_NumberStandardAATypes)
        );

        const size_t nr_angular_bins( ( *HISTOGRAM3DS)( 0)( 0)->GetNumberOfBinsY());
        const size_t nr_distance_bins( ( *HISTOGRAM3DS)( 0)( 0)->GetNumberOfBinsX());
        const double delta_x( ( *HISTOGRAM3DS)( 0)( 0)->GetBinSizeXYZ().First());
        const double start_x( ( *HISTOGRAM3DS)( 0)( 0)->GetBoundariesX().First());
        for( size_t type_a( 0), nr_aa_types( biol::AATypes::s_NumberStandardAATypes); type_a < nr_aa_types; ++type_a)
        {
          for( size_t type_b( 0); type_b <= type_a; ++type_b)
          {
            histogram_vec( type_a)( type_b) = histogram_vec( type_b)( type_a) =
              util::ShPtr< linal::Matrix< double> >( new linal::Matrix< double>( nr_angular_bins, nr_angular_bins, start_x + delta_x * nr_distance_bins));
            linal::Matrix< double> &matrix( *histogram_vec( type_a)( type_b));
            const math::Histogram3D &hist3d( *( *HISTOGRAM3DS)( type_a)( type_b));
            const math::Tensor< double> &tensor3d( hist3d.GetHistogram());
            for( size_t angle_a( 0); angle_a < nr_angular_bins; ++angle_a)
            {
              for( size_t angle_b( 0); angle_b < nr_angular_bins; ++angle_b)
              {
                for( size_t x_bin( 0); x_bin < nr_distance_bins; ++x_bin)
                {
                  if( tensor3d( x_bin, angle_a, angle_b) > 0.0)
                  {
                    matrix( angle_a, angle_b) = start_x + delta_x * x_bin;
                    break;
                  }
                }
              }
            }
          }
        }
        return histogram_vec;
      }
    }

    //! @brief read map of amino acid pair energies based on distance from histogram files
    void AAPairHiResClash::ReadEnergyFunctionMap()
    {
      m_Histograms = util::ToSiPtr( GetHistograms( m_HistogramFileName));
      m_DistanceCutoff = 5.5;
      m_ClashDistances = util::ToSiPtr( GetClashDistances( m_Histograms));
      m_AngularBinSize = ( *m_Histograms)( 0)( 0)->GetBinSizeXYZ()( 1);
    }

    //! @brief set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @return ERROR_STREAM stream with which to write errors
    bool AAPairHiResClash::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      if( !command::CommandState::IsInStaticInitialization())
      {
        ReadEnergyFunctionMap();
      }
      return true;
    }

    //! @brief Get histograms from a particular file. Caches histograms so they need only be read in once
    const storage::Vector< util::ShPtrVector< math::Histogram3D> > &AAPairHiResClash::GetHistograms( const std::string &FILENAME)
    {
      static storage::Map< std::string, storage::Vector< util::ShPtrVector< math::Histogram3D> > > s_histograms;
      storage::Vector< util::ShPtrVector< math::Histogram3D> > &histogram_vec( s_histograms[ FILENAME]);
      if( !histogram_vec.IsEmpty() || command::CommandState::IsInStaticInitialization())
      {
        return histogram_vec;
      }

      // store map of all aa pairs with their histogram
      histogram_vec.Reset();
      histogram_vec.Resize
      (
        biol::AATypes::s_NumberStandardAATypes,
        util::ShPtrVector< math::Histogram3D>( biol::AATypes::s_NumberStandardAATypes)
      );

      // read file with all histograms for each pair of aa types
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( FILENAME));

      // initialize temporary strings to read
      std::string tmp_a, tmp_b;

      // while reading the aatype pair and until reaching the end of the file
      while( read >> tmp_a >> tmp_b && !read.eof())
      {
        // read the string and convert it to first aa type
        biol::AAType a( biol::GetAATypes().AATypeFromOneLetterCode( tmp_a[0]));
        biol::AAType b( biol::GetAATypes().AATypeFromOneLetterCode( tmp_b[0]));

        // read the histogram
        util::ShPtr< math::Histogram3D> histogram( new math::Histogram3D);
        read >> *histogram;
        histogram_vec( a)( b) = histogram_vec( b)( a) = histogram;
      }

      // close the stream
      io::File::CloseClearFStream( read);
      return histogram_vec;
    }

  } // namespace score
} // namespace bcl
