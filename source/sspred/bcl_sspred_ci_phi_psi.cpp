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
#include "sspred/bcl_sspred_ci_phi_psi.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_chain.h"
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_aa_classes.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "coord/bcl_coord_polygon.h"
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_const_graph.h"
#include "linal/bcl_linal_vector_2d_operations.h"
#include "math/bcl_math_limits.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_running_min_max.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "sspred/bcl_sspred_pdb.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {

    //! @brief Data as string
    //! @param DATA the data whose name is desired
    //! @return the name as string
    const std::string &CIPhiPsi::GetTMDirectionType( const TMDirection &DATA)
    {
      static const std::string s_directions[ s_NumberTMDirections + 1] =
      {
        ">",
        "<",
        "a",
        "p",
        "i",
        "o",
        "-",
        "X",
        GetStaticClassName< CIPhiPsi::TMDirection>()
      };
      return s_directions[ DATA];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CIPhiPsi::CIPhiPsi() :
      m_SSType(),
      m_EnvironmentType(),
      m_TMDirection( e_NonMembrane),
      m_BetaBarrelResidueFacesPore( false),
      m_SSENumber( util::GetUndefinedSize_t())
    {
    }

    //! @brief constructor from members
    //! @param SS_TYPE type of sse
    //! @param ENVIRONMENT environment of the residue
    //! @param FACES_PORE true if the residue faces the pore of the membrane
    CIPhiPsi::CIPhiPsi
    (
      const biol::SSType &SS_TYPE,
      const biol::EnvironmentType &ENVIRONMENT,
      const bool &FACES_PORE,
      const TMDirection &DIRECTION
    ) :
      m_SSType( SS_TYPE),
      m_EnvironmentType( ENVIRONMENT),
      m_TMDirection( DIRECTION),
      m_BetaBarrelResidueFacesPore( FACES_PORE),
      m_SSENumber( util::GetUndefinedSize_t())
    {
    }

    //! @brief Clone function
    //! @return pointer to new CIPhiPsi
    CIPhiPsi *CIPhiPsi::Clone() const
    {
      return new CIPhiPsi( *this);
    }

    CIPhiPsi::~CIPhiPsi()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CIPhiPsi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &CIPhiPsi::GetFileExtension() const
    {
      static const std::string s_file_extension( ".ciphipsi");

      // end
      return s_file_extension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get whether the residue faces the pore
    //! @return true if the residue is in a beta barrel and faces the pore
    bool CIPhiPsi::DoesBetaBarrelResidueFacePore() const
    {
      return m_BetaBarrelResidueFacesPore;
    }

    //! @brief get whether the residue is in a TM-region whose N-term is towards the cytosol
    //! @return true if the residue is in a TM-region whose N-term is towards the cytosol
    CIPhiPsi::TMDirectionTypeEnum CIPhiPsi::GetTMDirection() const
    {
      return m_TMDirection;
    }

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D CIPhiPsi::GetThreeStatePrediction() const
    {
      return m_SSType->GetThreeStatePrediction();
    }

    //! @brief get the SSE Number
    //! @return the sse number
    const size_t &CIPhiPsi::GetSSENumber() const
    {
      return m_SSENumber;
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> CIPhiPsi::GetNineStatePrediction() const
    {
      return ConvertThreeStateToNineState( GetThreeStatePrediction(), m_EnvironmentType);
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &CIPhiPsi::ReadPredictionsForAA( std::istream &ISTREAM, biol::AABase &AMINO_ACID) const
    {
      std::string temp;

      // read chain id
      char chain_id;
      ISTREAM >> chain_id;
      while( ISTREAM.good() && chain_id < AMINO_ACID.GetChainID())
      {
        // call getline to read the next line
        std::getline( ISTREAM, temp);
        ISTREAM >> chain_id;
      }
      BCL_Assert
      (
        chain_id == AMINO_ACID.GetChainID() || ( AMINO_ACID.GetChainID() == ' ' && chain_id == '_'),
        "Chain ID mismatch; CIPhiPsi file had " + util::Format()( chain_id)
        + " but AA had " + util::Format()( AMINO_ACID.GetChainID())
      );

      // read secondary structure
      char ss_char;
      biol::SSType ss;
      ISTREAM >> ss_char;
      switch( ss_char)
      {
        case 'H': ss = biol::GetSSTypes().HELIX; break;
        case 'E': ss = biol::GetSSTypes().STRAND; break;
        case 'C': ss = biol::GetSSTypes().COIL; break;
        case 'X': ss = biol::GetSSTypes().e_Undefined; break;
        default:
          BCL_Exit( "Unrecognized SS Type: " + util::Format()( ss_char), -1);
      }

      // read membrane topology
      char env_char;
      ISTREAM >> env_char;
      const biol::EnvironmentType env_type
      (
        env_char == 'M'
        ? biol::GetEnvironmentTypes().e_MembraneCore
        : biol::GetEnvironmentTypes().e_Solution
      );

      // read orientation
      char orientation_char;
      ISTREAM >> orientation_char;
      const bool faces_pore( orientation_char == 'P');

      // call getline to read the next line
      std::getline( ISTREAM, temp);
      temp = util::TrimString( temp);
      if( temp.empty())
      {
        BCL_MessageStd( "Warning: Old CiPhiPsi file used. This is only valid if this is not a membrane protein!");
      }
      TMDirectionTypeEnum direction( temp);

      // determine method type
      const Method method( GetMethods().e_CIPhiPsi);
      AMINO_ACID.SetSSPrediction( method, CIPhiPsi( ss, env_type, faces_pore, direction));

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &CIPhiPsi::ReadPredictionsForAASequence( std::istream &ISTREAM, biol::AASequence &AA_SEQUENCE) const
    {
      // remove header lines
      util::ChopHeader( ISTREAM);

      // determine method type
      const Method method( GetMethods().e_CIPhiPsi);

      // get the desired chain id
      const char desired_chain_id( AA_SEQUENCE.GetChainID() == ' ' ? '_' : AA_SEQUENCE.GetChainID());

      // read in the seqid and one letter code
      int seqid;
      std::string one_letter_code, temp;
      char actual_chain_id;
      ISTREAM >> seqid >> one_letter_code >> actual_chain_id;
      while( ISTREAM.good() && actual_chain_id != desired_chain_id)
      {
        // read the rest of the line; then the next chain id
        std::getline( ISTREAM, temp);
        ISTREAM >> seqid >> one_letter_code >> actual_chain_id;
      }
      BCL_Assert
      (
        actual_chain_id == desired_chain_id || AA_SEQUENCE.Sequence().find_first_not_of( "XU") == std::string::npos,
        "Requested chain " + util::Format()( desired_chain_id) + " was not in CIPhiPsi file with sequence: " + AA_SEQUENCE.Sequence()
      );
      if( actual_chain_id != desired_chain_id)
      {
        return ISTREAM;
      }

      // get current aa type from read one letter code
      const biol::AAType current_aa_type( biol::GetAATypes().AATypeFromOneLetterCode( one_letter_code[ 0]));

      // get the first aa out of the sequence
      biol::AABase &first_aa( **AA_SEQUENCE.Begin());

      // assert matching amino acid seqids
      BCL_Assert
      (
        first_aa.GetSeqID() == seqid,
        "mismatch in seqids!\n sequence: " + first_aa.GetIdentification() +
        " vs. from sspred CIPhiPsi: " + util::Format()( seqid) + " " + one_letter_code
      );

      // check matching amino acid types
      BCL_Assert
      (
        !current_aa_type->IsNaturalAminoAcid() || !first_aa.GetType()->IsNaturalAminoAcid() ||
        current_aa_type == first_aa.GetType(),
        "mismatch in amino acid types! sequence: "
        + first_aa.GetIdentification()
        + " vs. from sspred: " + util::Format()( seqid) + " " + one_letter_code
      );
      ISTREAM.putback( actual_chain_id);
      this->ReadPredictionsForAA( ISTREAM, first_aa);

      // iterate over all amino acids in the sequence
      for
      (
        biol::AASequence::iterator aa_itr( AA_SEQUENCE.Begin() + 1), aa_itr_end( AA_SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // read the predictions for this amino acid
        MethodHandler::ReadPredictionsForAA( ISTREAM, **aa_itr, method);
      }

      SetSSENumbers( AA_SEQUENCE);

      // call standard read function and return it
      return ISTREAM;
    }

    namespace
    {

      //! @brief create a hull, initially convex, but subdividing each edge longer than the specified distance so as to
      //!        reach a polygon where the maximum distance of any side is no longer than specified
      //! @note this is useful for creating the bounding polygon for beta barrels because it prevents kinks in the barrel
      //!       from masking nearby strands
      //! @param POINTS the initial set of points
      //! @param MAX_SIDE_DISTANCE the maximum distance to allow for a given side of the polygon
      //! @return a polygon with the desired properties, if one exists
      coord::Polygon FindHull
      (
        const storage::Vector< linal::Vector2D> &POINTS,
        const double &MAX_SIDE_DISTANCE
      )
      {
        const coord::Polygon initial_hull( coord::Polygon::ConvexHull( POINTS));

        // determine which lines are longer than the desired distance
        const size_t initial_number_lines( initial_hull.GetSides().GetSize());

        // hash vector to determine which lines are to be split
        std::string line_must_be_split( initial_number_lines, '0');

        // determine which lines must be split
        size_t number_lines_to_split( 0);
        for( size_t line_number( 0); line_number < initial_number_lines; ++line_number)
        {
          if( initial_hull.GetSides()( line_number).GetLength() > MAX_SIDE_DISTANCE)
          {
            line_must_be_split[ line_number] = '1';
            ++number_lines_to_split;
          }
        }

        // if no lines were longer, return the convex hull
        if( !number_lines_to_split || initial_number_lines == size_t( 1))
        {
          return initial_hull;
        }

        // for each point, determine which line it lies nearest
        // and place the point into a container for that line, if the line is to be split
        storage::Vector< storage::Vector< linal::Vector2D> > points_nearest_line( initial_number_lines);

        // get the barycenter of the initial hull
        const linal::Vector2D initial_barycenter( initial_hull.GetBarycenter());
        for
        (
          storage::Vector< linal::Vector2D>::const_iterator itr_points( POINTS.Begin()), itr_points_end( POINTS.End());
          itr_points != itr_points_end;
          ++itr_points
        )
        {
          // skip corner points, since they are nearest to two different lines and must be the start and end points of
          // any new paths
          if( initial_hull.IsCornerOf( *itr_points))
          {
            continue;
          }

          // find the closest line and distance from it
          std::pair< coord::Polygon::const_iterator, double> closest_line_and_distance
          (
            initial_hull.FindNearestSide( *itr_points)
          );

          const size_t nearest_line_id( std::distance( initial_hull.Begin(), closest_line_and_distance.first));
          if( line_must_be_split[ nearest_line_id] == '1')
          {
            points_nearest_line( nearest_line_id).PushBack( *itr_points);
          }
        }

        // create a new polygon to contain the distance limited hull
        coord::Polygon distance_limited_hull;

        // for each line
        const double undefined_edge_value( math::GetHighestBoundedValue< double>());
        for( size_t line_number( 0); line_number < initial_number_lines; ++line_number)
        {
          const coord::LineSegment2D &current_side( initial_hull.GetSides()( line_number));
          if( line_must_be_split[ line_number] == '0' || points_nearest_line( line_number).IsEmpty())
          {
            distance_limited_hull.PushBack( current_side.GetStartPoint());
            continue;
          }
          // add the starting and ending points to the nearest lines
          points_nearest_line( line_number).PushBack( current_side.GetStartPoint());
          points_nearest_line( line_number).PushBack( current_side.GetEndPoint());
          const storage::Vector< linal::Vector2D> &nearest_points( points_nearest_line( line_number));

          // determine the indices of the two endpoints on this hull line in the points_nearest_line vector
          const size_t start_point_id( nearest_points.GetSize() - 2);
          const size_t end_point_id( nearest_points.GetSize() - 1);

          // get the path with small side lengths that goes from one side of the convex hull on this side to the other
          // this is seeking to minimize the sum of areas given by rectangles drawn around each line
          graph::Path path;
          const size_t n_points( nearest_points.GetSize());
          graph::ConstGraph< size_t, double> distance_graph;
          {
            linal::Matrix< double> distances( n_points, n_points, undefined_edge_value);
            for( size_t point_a( 0); point_a < n_points; ++point_a)
            {
              for( size_t point_b( point_a + 1); point_b < n_points; ++point_b)
              {
                distances( point_a, point_b)
                  = linal::SquareDistance( nearest_points( point_a), nearest_points( point_b));
              }
            }

            distance_graph = graph::ConstGraph< size_t, double>
                             (
                               storage::Vector< size_t>( n_points, size_t( 1)),
                               distances,
                               undefined_edge_value,
                               false,
                               true
                             );
          }
          path = graph::Connectivity::FindMinimalPath( distance_graph, start_point_id, end_point_id);

          if( path.GetSize() <= size_t( 2))
          {
            distance_limited_hull.PushBack( current_side.GetStartPoint());
            continue;
          }

          // create a polygon from the given path, skipping the last point which is the start point for the next path
          for( graph::Path::const_iterator itr( path.Begin()), itr_end( path.End() - 1); itr != itr_end; ++itr)
          {
            distance_limited_hull.PushBack( nearest_points( *itr));
          }
        }
        return distance_limited_hull;
      } // end CreateConvexHullWithLimitedDistanceSides
    }

    //! @brief iterates over the sequences in ProteinModel and calculates CIPhiPsi for every residue in the sequence
    //! @param PROTEIN_MODEL ProteinModel for which CIPhiPsi will be calculated
    void CIPhiPsi::Calculate( assemble::ProteinModel &PROTEIN_MODEL, const bool &USE_PDBTM_MEMBRANE_THICKNESS) const
    {
      // get membrane
      util::ShPtr< biol::Membrane> sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      // if membrane is not defined, try to read it in
      if( !sp_membrane.IsDefined())
      {
        PDB::SetEnvironmentTypes( PROTEIN_MODEL, USE_PDBTM_MEMBRANE_THICKNESS);
        sp_membrane = PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane);
      }

      // iterate through the chains
      for
      (
        util::ShPtrVector< assemble::Chain>::iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        biol::AASequence &sequence( *( *chain_itr)->GetSequence());
        // call calculate on the sequence
        Calculate( sequence, util::SiPtr< const biol::Membrane>( sp_membrane));
      }
    }

    //! @brief iterates over the sequence and calculates CIPhiPsi for every residue in the sequence
    //! @param SEQUENCE Sequence of interest
    //! @param MEMBRANE_PTR pointer to membrane object
    void CIPhiPsi::Calculate( biol::AASequence &SEQUENCE, const util::SiPtr< const biol::Membrane> &MEMBRANE_PTR) const
    {
      // set the secondary structure according to basic angular information
      CalculateSolubleSS
      (
        SEQUENCE,
        MEMBRANE_PTR.IsDefined() ? MEMBRANE_PTR->GetThickness( biol::GetEnvironmentTypes().e_MembraneCore) : 0.0
      );

      if( !MEMBRANE_PTR.IsDefined())
      {
        return;
      }

      // determine method type
      const Method method( GetMethods().e_CIPhiPsi);

      // Only update residues in the membrane core region, so determine how thick that region is
      const biol::Membrane &membrane( *MEMBRANE_PTR);
      const double membrane_thickness( membrane.GetLimit( biol::GetEnvironmentTypes().e_MembraneCore));
      if( membrane_thickness == 0.0)
      {
        // no membrane core -> soluble protein
        return;
      }

      // slice the membrane into approximately 3A slices, store coordinates of each point
      const double approximate_slice_thickness( 3.8);

      // adjust the slice thickness such that it is as near as possible to the slice thickness
      // while still being an integral multiple of 2A
      const double slice_thickness
      (
        2.0 * membrane_thickness
        /
        double( std::max( int( membrane_thickness / ( approximate_slice_thickness / 2.0)), int( 1)))
      );

      // number of slices
      const size_t number_slices( 2.0 * membrane_thickness / slice_thickness);

      // create a vector with the aa's in each slice
      storage::Vector< storage::List< util::SiPtr< biol::AABase> > > aas_in_slice( number_slices);

      // create a list of all membrane helix residues
      storage::List< util::SiPtr< biol::AABase> > membrane_helix_residues;

      // raw coordinates for each slice
      storage::Vector< storage::Vector< linal::Vector2D> > slice_coordinates( number_slices);

      // iterate through the sequence
      size_t membrane_helix_residue_counts( 0), membrane_strand_residue_counts( 0), membrane_coil_residue_counts( 0);
      biol::SSType previous_ss( biol::GetSSTypes().COIL);
      int previous_ss_pdb_id_start( math::GetLowestBoundedValue< int>());
      double previous_z_coord_start( util::GetUndefined< double>());
      size_t tm_strand_count( 0), tm_helix_count( 0);

      // iterate over all residues in the given sequence
      for
      (
        biol::AASequence::iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // get a const-ref to the coordinates of this residues' CA
        const linal::Vector3D &coordinates( ( *aa_itr)->GetCA().GetCoordinates());

        // determine the z-coordinate of this aa
        const double z_coordinate( coordinates.Z());
        if( !util::IsDefined( z_coordinate))
        {
          continue;
        }

        // test whether this residue has a ss identification
        util::SiPtr< const CIPhiPsi> ciphipsi_sol( ( *aa_itr)->GetSSPrediction( method));
        if( !ciphipsi_sol.IsDefined())
        {
          continue;
        }

        // ciphipsi only differs from stridedssp in the membrane core region, so skip residues that do not appear to be
        // in the core
        const double abs_z_coordinate( math::Absolute( z_coordinate));

        // if the residue is clearly in the solution, continue
        if( abs_z_coordinate >= membrane_thickness)
        {
          if( previous_ss == biol::GetSSTypes().STRAND)
          {
            if( ( *aa_itr)->GetPdbID() - previous_ss_pdb_id_start >= 6)
            {
              ++tm_strand_count;
            }
          }
          else if( previous_ss == biol::GetSSTypes().HELIX)
          {
            if( ( *aa_itr)->GetPdbID() - previous_ss_pdb_id_start >= 10)
            {
              ++tm_helix_count;
            }
          }
          previous_ss = biol::GetSSTypes().COIL;
          continue;
        }

        // get ss prediction
        const biol::SSType ss( ciphipsi_sol->GetOneStateSSPrediction());

        if( previous_ss != ss)
        {
          if( previous_ss == biol::GetSSTypes().STRAND)
          {
            if( ( *aa_itr)->GetPdbID() - previous_ss_pdb_id_start >= 6)
            {
              ++tm_strand_count;
            }
          }
          else if( previous_ss == biol::GetSSTypes().HELIX)
          {
            if( ( *aa_itr)->GetPdbID() - previous_ss_pdb_id_start >= 10)
            {
              ++tm_helix_count;
            }
          }
          if( ss != biol::GetSSTypes().COIL)
          {
            previous_ss_pdb_id_start = ( *aa_itr)->GetPdbID();
          }
          previous_ss = ss;
        }

        if( ss == biol::GetSSTypes().HELIX)
        {
          membrane_helix_residues.PushBack( **aa_itr);
        }

        // update the CIPhiPsi analysis to denote that this is a membrane core residue
        // This will only be changed if the residue is detected to be in a pore
        if( ss == biol::GetSSTypes().COIL)
        {
          ( *aa_itr)->SetSSPrediction
          (
            method,
            CIPhiPsi( ss, biol::GetEnvironmentTypes().e_MembraneCore, false, ciphipsi_sol->GetTMDirection())
          );
        }
        else
        {
          ( *aa_itr)->SetSSPrediction
          (
            method,
            CIPhiPsi( ss, biol::GetEnvironmentTypes().e_MembraneCore, false, ciphipsi_sol->GetTMDirection())
          );
        }

        // determine which slice id this residue belongs in
        const size_t slice_id
        (
          std::min( size_t( ( z_coordinate + membrane_thickness) / slice_thickness), number_slices - 1)
        );

        if( ss == biol::GetSSTypes().COIL)
        {
          ++membrane_coil_residue_counts;
        }
        else if( ss == biol::GetSSTypes().STRAND)
        {
          ++membrane_strand_residue_counts;
        }
        else if( ss == biol::GetSSTypes().HELIX)
        {
          ++membrane_helix_residue_counts;
          continue;
        }

        //BCL_MessageStd( ( *aa_itr)->GetIdentification() + " slice id: " + util::Format()( slice_id) + " " + util::Format()( abs_z_coordinate));

        // append this residue to that slice
        slice_coordinates( slice_id).PushBack( linal::Vector2D( coordinates.X(), coordinates.Y()));
        aas_in_slice( slice_id).PushBack( **aa_itr);
      }
      if( previous_ss == biol::GetSSTypes().STRAND)
      {
        if( ( *--SEQUENCE.End())->GetPdbID() - previous_ss_pdb_id_start >= 6)
        {
          ++tm_strand_count;
        }
      }
      else if( previous_ss == biol::GetSSTypes().HELIX)
      {
        if( ( *--SEQUENCE.End())->GetPdbID() - previous_ss_pdb_id_start >= 10)
        {
          ++tm_helix_count;
        }
      }

      // determine the total number of residues in the membrane
      const size_t membrane_residue_counts
      (
        membrane_helix_residue_counts + membrane_strand_residue_counts + membrane_coil_residue_counts
      );

      // detect whether this really appears to be a beta barrel protein
      if( tm_helix_count > tm_strand_count || membrane_residue_counts < 16)
      {
        // not a beta barrel, so no pore.  Thus there is no need to refine analysis by looking for BB pores
        return;
      }

      // is a beta barrel.  Change all membrane helices to soluble helices
      for
      (
        storage::List< util::SiPtr< biol::AABase> >::iterator
          aa_itr( membrane_helix_residues.Begin()), aa_itr_end( membrane_helix_residues.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // set the helix to be in the solution; likewise, revert to StrideDSSP prediction (usually also helix)
        ( *aa_itr)->SetSSPrediction
        (
          method,
          CIPhiPsi
          (
            ( *aa_itr)->GetSSPrediction( method)->GetOneStateSSPrediction(),
            biol::GetEnvironmentTypes().e_Solution,
            false,
            e_Pore
          )
        );
      }

      BCL_MessageVrb
      (
        "Membrane strand residues: " + util::Format()( membrane_strand_residue_counts)
        + " # tm strands: " + util::Format()( tm_strand_count)
        + " # elements: " + util::Format()( SEQUENCE.GetSize())
      );

      // track all coordinates in the membrane
      storage::Vector< linal::Vector2D> membrane_coordinates;

      // combine adjacent slices such that each slice has a minimum of one residue per strand, on average
      for
      (
        size_t slice_id( 0), number_slices_before_end( number_slices - 1);
        slice_id < number_slices_before_end;
        ++slice_id
      )
      {
        if( slice_coordinates( slice_id).GetSize() < tm_strand_count)
        {
          // combine with the next slice, if possible
          slice_coordinates( slice_id + 1).Append( slice_coordinates( slice_id));
          slice_coordinates( slice_id).Reset();
          aas_in_slice( slice_id + 1).Append( aas_in_slice( slice_id));
          aas_in_slice( slice_id).Reset();
        }
        else
        {
          membrane_coordinates.Append( slice_coordinates( slice_id));
        }
      }
      membrane_coordinates.Append( slice_coordinates.LastElement());

      // if the last slice has anything left in it < tm_strand_count, combine with previous slices
      for( size_t slice_id( number_slices - 1); slice_id; --slice_id)
      {
        if( slice_coordinates( slice_id).GetSize() < tm_strand_count)
        {
          // combine with the next slice, if possible
          slice_coordinates( slice_id - 1).Append( slice_coordinates( slice_id));
          slice_coordinates( slice_id).Reset();
          aas_in_slice( slice_id - 1).Append( aas_in_slice( slice_id));
          aas_in_slice( slice_id).Reset();
        }
      }

      // compute the barycenter of the beta barrel
      const linal::Vector2D membrane_barycenter( coord::Polygon::ConvexHull( membrane_coordinates).GetBarycenter());

      // for each slice, find all residues within 3A of the convex hull
      for( size_t slice_id( 0); slice_id < number_slices; ++slice_id)
      {
        // get the coordinates for this slice
        const storage::Vector< linal::Vector2D> &this_slice_coordinates( slice_coordinates( slice_id));

        // skip slices that are empty
        if( this_slice_coordinates.IsEmpty())
        {
          continue;
        }

        // get the high-res hull
        const coord::Polygon hull( FindHull( this_slice_coordinates, 7.0));

        // declare all AAs within 3 angstroms of a side of the hull to be in the membrane; others are in a pore
        // for those in the membrane, determine orientation
        storage::List< util::SiPtr< biol::AABase> >::iterator
          itr_aa( aas_in_slice( slice_id).Begin());
        for
        (
          storage::Vector< linal::Vector2D>::const_iterator
            itr_slice_coords( this_slice_coordinates.Begin()), itr_slice_coords_end( this_slice_coordinates.End());
          itr_slice_coords != itr_slice_coords_end;
          ++itr_slice_coords, ++itr_aa
        )
        {
          if( hull.FindNearestSide( *itr_slice_coords).second > 3.8)
          {
            // point in pore; reset to StrideDSSP prediction
            ( *itr_aa)->SetSSPrediction
            (
              method,
              CIPhiPsi
              (
                ( *itr_aa)->GetSSPrediction( method)->GetOneStateSSPrediction(),
                biol::GetEnvironmentTypes().e_Solution,
                false,
                e_Pore
              )
            );
            continue;
          }

          // aa is in the membrane, update with orientation information, if it is available

          const linal::Vector3D &first_sidechain_3d( ( *itr_aa)->GetFirstSidechainAtom().GetCoordinates());
          linal::Vector2D first_sidechain_2d( first_sidechain_3d.X(), first_sidechain_3d.Y());
          math::RunningAverage< linal::Vector2D> average_sidechain_position;
          average_sidechain_position += first_sidechain_2d;
          for
          (
            util::SiPtrVector< const biol::Atom>::const_iterator
              itr_atom( ( *itr_aa)->GetAtoms().Begin()), itr_atom_end( ( *itr_aa)->GetAtoms().End());
            itr_atom != itr_atom_end;
            ++itr_atom
          )
          {
            linal::Vector2D position_2d( ( *itr_atom)->GetCoordinates().X(), ( *itr_atom)->GetCoordinates().Y());
            if( ( *itr_atom)->GetType()->IsSideChain())
            {
              average_sidechain_position += position_2d;
            }
          }
          const bool is_closer
          (
            linal::Distance( average_sidechain_position.GetAverage(), membrane_barycenter)
            < linal::Distance( *itr_slice_coords, membrane_barycenter)
          );
          util::SiPtr< const CIPhiPsi> ciphipsi_si( ( *itr_aa)->GetSSPrediction( method));
          if( is_closer)
          {
            ( *itr_aa)->SetSSPrediction
            (
              method,
              CIPhiPsi
              (
                biol::GetSSTypes().STRAND,
                biol::GetEnvironmentTypes().e_MembraneCore,
                true,
                ciphipsi_si->GetTMDirection()
              )
            );
          }
          else
          {
            ( *itr_aa)->SetSSPrediction
            (
              method,
              CIPhiPsi
              (
                biol::GetSSTypes().STRAND,
                biol::GetEnvironmentTypes().e_MembraneCore,
                false,
                ciphipsi_si->GetTMDirection()
              )
            );
          }
        }
      }
      SetSSENumbers( SEQUENCE);
    }

    //! @brief iterates over the sequence and writes CIPhiPsi for every residue in every chain
    //! @param STREAM stream to write the analysis to
    //! @param PROTEIN_MODEL Model to write the analysis for
    void CIPhiPsi::WriteAnalysis( std::ostream &STREAM, const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // write the normal header
      STREAM << "# Context-insensitive secondary structure analysis file\n"
       "# Version 1.0\n"
       "# SecondaryStructure is one of H/E/C ( Helix / Strand / Coil )\n"
       "# MembraneLocalization is one of M/S ( Membrane / Solution )\n"
       "# Orientation is either L (residue faces lipid), P (residue faces pore), or -\n"
       "# Orientation of '-' is used for soluble residues\n"
       "# Origin is c for TM-SS elements originating in the cytosol, or e for extracellular/outside cytosol\n"
       "# SeqID 1LtrCode Chain SecondaryStructure MembraneLocalization Orientation Origin\n";

      // iterate through the chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // call WriteAnalysis on the sequence
        WriteAnalysis( STREAM, *( *chain_itr)->GetSequence(), ( *chain_itr)->GetChainID());
      }
    }

    //! @brief iterates over the protein and calculates CIPhiPsi for every residue in every chain
    //! @param STREAM stream to write the analysis to
    //! @param SEQUENCE Sequence to write the analysis for
    //! @param CHAIN_ID Id of the chain that the sequence represents
    void CIPhiPsi::WriteAnalysis( std::ostream &STREAM, const biol::AASequence &SEQUENCE, const char &CHAIN_ID) const
    {
      // determine method type
      const Method method( GetMethods().e_CIPhiPsi);
      const char chain_id( CHAIN_ID == ' ' ? '_' : CHAIN_ID);

      // iterate over amino acids in the sequence
      for
      (
        biol::AASequence::const_iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        util::SiPtr< const CIPhiPsi> ciphipsi_ptr( ( *aa_itr)->GetSSPrediction( method));

        // check existence of requested method
        if( !ciphipsi_ptr.IsDefined())
        {
          BCL_MessageCrt( "CIPhiPsi was not stored for " + ( *aa_itr)->GetIdentification());
          if( util::IsDefined( ( *aa_itr)->GetSeqID()))
          {
            STREAM << ( *aa_itr)->GetSeqID() << ' '
                   << ( *aa_itr)->GetType()->GetOneLetterCode() << ' '
                   << chain_id << " X S X X\n";
          }
          continue;
        }

        const CIPhiPsi &ciphipsi( *ciphipsi_ptr);

        // write residue info
        STREAM << ( *aa_itr)->GetSeqID() << ' ' << ( *aa_itr)->GetType()->GetOneLetterCode() << ' ' << chain_id << ' ';

        // write analysis info
        STREAM << ciphipsi.GetOneStateSSPrediction()->GetOneLetterCode() << ' ';
        if( ciphipsi.GetOneStateTMPrediction() == biol::GetEnvironmentTypes().e_MembraneCore)
        {
          STREAM << "M ";
          STREAM << ( ciphipsi.DoesBetaBarrelResidueFacePore() ? 'P' : 'L');
        }
        else
        {
          // soluble residue
          STREAM << "S -";
        }
        STREAM << ' ' << GetTMDirectionType( ciphipsi.GetTMDirection()) << '\n';
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CIPhiPsi::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSType, ISTREAM);
      io::Serialize::Read( m_EnvironmentType, ISTREAM);
      io::Serialize::Read( m_BetaBarrelResidueFacesPore, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CIPhiPsi::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_EnvironmentType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BetaBarrelResidueFacesPore, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

    //! @brief assign SS for all residues, ignoring the membrane
    //! @param SEQUENCE Sequence of interest
    //! @param MEMBRANE_PTR pointer to membrane object
    void CIPhiPsi::CalculateSolubleSS( biol::AASequence &SEQUENCE, const double &MEM_THICKNESS) const
    {
      // 1-3 residues cannot be categorized
      if( SEQUENCE.GetSize() < size_t( 4))
      {
        return;
      }

      // determine method type
      const Method method( GetMethods().e_CIPhiPsi);

      // create an assemble protein model with cache from the sequence
      assemble::ProteinModelWithCache pmwc
      (
        assemble::ProteinModel
        (
          util::ShPtr< assemble::Chain>
          (
            new assemble::Chain
            (
              util::ShPtr< biol::AASequence>( SEQUENCE.Clone())
            )
          )
        ),
        true
      );

      if( pmwc.GetSize() < size_t( 4))
      {
        return;
      }

      // this descriptor is really a decision tree that was trained to distinguish strands from coils
      // based on various phi-psi descriptors.
      // The interpretation is that a strand requires that there be at least one window of 3 residues
      // around an amino acid for which the maximum deviation from perfect strand phi-psi angles is less than
      // 50.59 degrees for both phi and psi for all residues within the window.
      // Further, a strand requires that either there be at least one window of 5 residues
      // around an amino acid for which the sum of deviations from perfect strand phi-psi angles is less than
      // 68.961 degrees for for all residues within the window OR the distance spanned by the nine-residue window
      // centered on this AA is at least 18.166 angstroms.
      // While this descriptor could be written out in pure logical form; it is truly a pain to get it all right
      util::Implementation< descriptor::Base< biol::AABase, float> > strand_determiner
      (
        "Multiply("
        "  LessEqual("
        "    rhs=Constant(50.5862),"
        "    lhs=DescriptorMin("
        "          ReflectingWindow("
        "            DescriptorMax("
        "              ReflectingWindow("
        "                DescriptorMax("
        "                  Combine("
        "                    Abs(Subtract(lhs=Psi(origin=-90),rhs=Constant(135))),"
        "                    Abs(Subtract(lhs=Phi(origin=-10),rhs=Constant(230)))"
        "                  )"
        "                ),"
        "                size=1,"
        "                alignment=Center"
        "              )"
        "            ),"
        "            size=1,"
        "            alignment=Center"
        "          )"
        "        )"
        "  ),"
        "  Add("
        "    LessEqual("
        "      rhs=Constant(68.961),"
        "      lhs=DescriptorMin("
        "            ReflectingWindow("
        "              DescriptorMax("
        "                ReflectingWindow("
        "                  Add("
        "                    Abs(Subtract(lhs=Psi(origin=-90),rhs=Constant(135))),"
        "                    Abs(Subtract(lhs=Phi(origin=-10),rhs=Constant(230)))"
        "                  ),"
        "                  size=2,"
        "                  alignment=Center"
        "                )"
        "              ),"
        "              size=2,"
        "              alignment=Center"
        "            )"
        "          )"
        "    ),"
        "    GreaterEqual("
        "      rhs=Constant(18.166),"
        "      lhs=Sqrt(DescriptorSum(Sqr(Subtract(lhs=Offset(Position,offset=-4),rhs=Offset(Position,offset=4)))))"
        "    )"
        "  )"
        ")"
      );

      strand_determiner->SetObject( pmwc);

      // find all stretches of residues in the sequence for which phi, psi, and CA coordinates are known
      // if any stretch is fewer than 4 in length, it is discarded
      // AA, phi, psi for all AAs in known_phi_psi_ca
      storage::Vector< storage::Vector< t_SegmentStorage> > phipsica;

      // store whether we are currently in a phi-psi-ca containing segment
      bool in_phi_psi_ca_segment( false);

      // get radians to degrees conversion
      const double degrees_per_radian( 180.0 / math::g_Pi);

      descriptor::Iterator< biol::AABase> itr_descriptor( pmwc.GetIterator());
      // iterate over all residues in the given sequence
      for
      (
        biol::AASequence::iterator
          aa_itr_prev( SEQUENCE.Begin()),
          aa_itr( ++SEQUENCE.Begin()),
          aa_itr_next( ++++SEQUENCE.Begin()),
          aa_itr_end( SEQUENCE.End());
        aa_itr_next != aa_itr_end;
        ++aa_itr, ++aa_itr_prev, ++aa_itr_next
      )
      {
        biol::AABase &aa( **aa_itr);

        const biol::AABase &aa_prev( **aa_itr_prev);
        const biol::AABase &aa_next( **aa_itr_next);

        const int aa_seq_id_next( aa_next.GetSeqID());
        const int aa_seq_id_prev( aa_prev.GetSeqID());

        if( aa_seq_id_next - aa_seq_id_prev != 2)
        {
          if( in_phi_psi_ca_segment && phipsica.GetSize() && phipsica.LastElement().GetSize() < size_t( 4))
          {
            phipsica.PopBack();
          }
          in_phi_psi_ca_segment = false;
          continue;
        }
        while( itr_descriptor.NotAtEnd() && itr_descriptor( 0)->GetSeqID() != aa.GetSeqID())
        {
          ++itr_descriptor;
        }
        if( !itr_descriptor.NotAtEnd())
        {
          itr_descriptor.GotoPosition( 0);
        }

        // compute phi/psi
        double phi( aa.CalculatePhi( aa_prev.GetAtom( biol::GetAtomTypes().C)) * degrees_per_radian);
        double psi( aa.CalculatePsi( aa_next.GetAtom( biol::GetAtomTypes().N)) * degrees_per_radian);

        // get a reference to the coordinates for CA
        const linal::Vector3D &ca_coords( aa.GetCA().GetCoordinates());

        // check that everything is defined
        if( !util::IsDefined( phi) || !util::IsDefined( psi) || !ca_coords.IsDefined())
        {
          if( in_phi_psi_ca_segment && phipsica.GetSize() && phipsica.LastElement().GetSize() < size_t( 4))
          {
            phipsica.PopBack();
          }
          in_phi_psi_ca_segment = false;
          continue;
        }

        // wrap phi around at -10 degrees and psi around at -90; this keeps the distributions of points on the
        // ramachandran plot relatively continuous
        if( phi < -10.0)
        {
          phi += 360.0;
        }
        if( psi < -90.0)
        {
          psi += 360.0;
        }

        if( !in_phi_psi_ca_segment)
        {
          in_phi_psi_ca_segment = true;
          phipsica.PushBack();
        }

        phipsica.LastElement().PushBack
        (
          t_SegmentStorage
          (
            util::SiPtr< biol::AABase>( aa),
            storage::VectorND< 2, double>( phi, psi),
            strand_determiner->operator()( itr_descriptor)( 0)
          )
        );
      }
      if( phipsica.GetSize() && phipsica.LastElement().GetSize() < size_t( 4))
      {
        phipsica.PopBack();
      }

      // iterate over all segments
      for
      (
        storage::Vector< storage::Vector< t_SegmentStorage> >::iterator
          itr( phipsica.Begin()), itr_end( phipsica.End());
        itr != itr_end;
        ++itr
      )
      {
        storage::Vector< t_SegmentStorage> &segment( *itr);
        const storage::Vector< biol::SSType> ss_types( CalculateSolubleSSForSegment( segment, MEM_THICKNESS));
        const storage::Vector< TMDirectionTypeEnum> tm_direction_types
        (
          CalculateTMDirectionsForSegment( segment, MEM_THICKNESS)
        );
        storage::Vector< biol::SSType>::const_iterator itr_ss( ss_types.Begin());
        auto itr_tm( tm_direction_types.Begin());
        for
        (
          storage::Vector< t_SegmentStorage>::iterator
            itr_segment( segment.Begin()), itr_segment_end( segment.End());
          itr_segment != itr_segment_end;
          ++itr_segment, ++itr_ss, ++itr_tm
        )
        {
          itr_segment->First()->SetSSPrediction
          (
            method,
            CIPhiPsi
            (
              *itr_ss,
              *itr_tm == e_NonMembrane
              ? biol::GetEnvironmentTypes().e_Solution
              : biol::GetEnvironmentTypes().e_MembraneCore,
              false,
              *itr_tm
            )
          );
        }
      }

      // now iterate over and try to fill in gaps
      bool updated_previous( false);
      CIPhiPsi undefined_ciphipsi( biol::GetSSTypes().e_Undefined, biol::GetEnvironmentTypes().e_Undefined, false, e_Unknown);
      for
      (
        biol::AASequence::iterator
          aa_itr_prev( SEQUENCE.Begin()),
          aa_itr( ++SEQUENCE.Begin()),
          aa_itr_next( ++++SEQUENCE.Begin()),
          aa_itr_first( SEQUENCE.Begin()),
          aa_itr_last( --SEQUENCE.End()),
          aa_itr_end( SEQUENCE.End());
        aa_itr_next != aa_itr_end;
        ++aa_itr, ++aa_itr_prev, ++aa_itr_next
      )
      {
        util::SiPtr< const MethodInterface> pred( ( *aa_itr)->GetSSPrediction( method));
        int seq_id_aa_itr( ( *aa_itr)->GetSeqID());
        int seq_id_aa_itr_prev( ( *aa_itr_prev)->GetSeqID());
        int seq_id_aa_itr_next( ( *aa_itr_next)->GetSeqID());

        if( pred.IsDefined())
        {
          // handle edge cases (first and last aa)
          if( aa_itr_prev == aa_itr_first && seq_id_aa_itr == seq_id_aa_itr_prev + 1)
          {
            ( *aa_itr_prev)->SetSSPrediction
            (
              method,
              *util::SiPtr< const CIPhiPsi>( pred)
            );
          }
          if( aa_itr_next == aa_itr_last && seq_id_aa_itr == seq_id_aa_itr_next - 1)
          {
            ( *aa_itr_next)->SetSSPrediction
            (
              method,
              *util::SiPtr< const CIPhiPsi>( pred)
            );
          }
          continue;
        }
        if( updated_previous)
        {
          updated_previous = false;
          continue;
        }
        util::SiPtr< const MethodInterface> pred_prev
        (
          seq_id_aa_itr == seq_id_aa_itr_prev + 1
          ? ( *aa_itr_prev)->GetSSPrediction( method)
          : util::SiPtr< const MethodInterface>()
        );
        util::SiPtr< const MethodInterface> pred_next
        (
          seq_id_aa_itr == seq_id_aa_itr_next - 1
          ? ( *aa_itr_next)->GetSSPrediction( method)
          : util::SiPtr< const MethodInterface>()
        );
        // if both are defined or undefined, nothing left to do
        if( !pred_prev.IsDefined() && !pred_next.IsDefined())
        {
          continue;
        }
        util::SiPtr< const CIPhiPsi> pred_derived( pred_prev.IsDefined() ? pred_prev : pred_next);
        ( *aa_itr)->SetSSPrediction( method, *pred_derived);
        updated_previous = true;
      }

      SetSSENumbers( SEQUENCE);
    }

    //! @brief set SSE Numbers for a given AA Sequence
    void CIPhiPsi::SetSSENumbers( biol::AASequence &AA_SEQUENCE)
    {
      size_t current_sse_number( 1);
      const Method ciphipsi( GetMethods().e_CIPhiPsi);
      int prev_seq_id( AA_SEQUENCE.GetFirstAA()->GetSeqID() - 1);
      biol::SSType prev_prediction( biol::GetSSTypes().COIL);
      biol::SSType current_prediction( biol::GetSSTypes().COIL);
      bool was_in_membrane( false);
      std::string going_up_str, length_str;

      for
      (
        biol::AASequence::iterator itr_seq( AA_SEQUENCE.Begin()), itr_seq_end( AA_SEQUENCE.End());
        itr_seq != itr_seq_end;
        ++itr_seq
      )
      {
        util::SiPtr< const CIPhiPsi> sstype( ( *itr_seq)->GetSSPrediction( ciphipsi));
        if( !sstype.IsDefined())
        {
          continue;
        }
        current_prediction = sstype->GetOneStateSSPrediction();
        if( ( *itr_seq)->GetSeqID() != prev_seq_id + 1 || prev_prediction != current_prediction)
        {
          ++current_sse_number;
        }
        auto tm_type( sstype->GetOneStateTMPrediction());
        bool in_membrane
        (
          (
            tm_type == biol::GetEnvironmentTypes().e_MembraneCore
            || tm_type == biol::GetEnvironmentTypes().e_MembraneInside
            || tm_type == biol::GetEnvironmentTypes().e_MembraneOutside
          )
          && current_prediction->IsStructured()
        );
        CIPhiPsi newsstype( *sstype);
        newsstype.m_SSENumber = current_sse_number;
        if( in_membrane && !was_in_membrane)
        {
          // look ahead to the end of the membrane segment. If it is above this point, set going up flag
          biol::AASequence::iterator itr_seq_b( itr_seq);
          for( ; itr_seq_b != itr_seq_end; ++itr_seq_b)
          {
            util::SiPtr< const CIPhiPsi> sstype_b( ( *itr_seq_b)->GetSSPrediction( ciphipsi));
            if( !sstype_b.IsDefined())
            {
              continue;
            }
            if
            (
              sstype_b->GetOneStateTMPrediction() != tm_type
              || (
                   sstype_b->GetOneStateSSPrediction() != current_prediction
                   && ( itr_seq_b - itr_seq) >= int( 3)
                 )
            )
            {
              --itr_seq_b;
              break;
            }
          }
          if( itr_seq_b == itr_seq_end)
          {
            --itr_seq_b;
          }
          if( itr_seq_b == itr_seq)
          {
            // lone residue in the membrane. Delete
            in_membrane = false;
            newsstype.m_EnvironmentType = biol::GetEnvironmentTypes().e_Solution;
          }
          else
          {
            if( current_prediction == biol::GetSSTypes().HELIX)
            {
              BCL_MessageStd
              (
                "Debug #Res: " + util::Format()( size_t( itr_seq_b - itr_seq))
                + " " + util::Format()( ( *itr_seq)->GetSeqID()) + "-" + util::Format()( ( *itr_seq_b)->GetSeqID())
                + " Z:" + util::Format()( ( ( *itr_seq)->GetCA().GetCoordinates().Z())) + " - " + util::Format()( ( ( *itr_seq_b)->GetCA().GetCoordinates().Z()))
                + " " + util::Format()( ( ( *itr_seq_b)->GetCA().GetCoordinates().Z() - ( *itr_seq)->GetCA().GetCoordinates().Z()))
                + " " + util::Format()( ( ( *itr_seq_b)->GetCA().GetCoordinates().Z() - ( *itr_seq)->GetCA().GetCoordinates().Z()) / ( float( itr_seq_b - itr_seq) + 1.0))
              );
            }
            going_up_str += newsstype.GetTMDirection().GetString() + "  ";
            length_str += util::Format().W( 3)( size_t( itr_seq_b - itr_seq) + 1);
          }
        }
        was_in_membrane = in_membrane;

        ( *itr_seq)->SetSSPrediction( ciphipsi, newsstype);

        util::SiPtr< const CIPhiPsi> sstypen( ( *itr_seq)->GetSSPrediction( ciphipsi));
        prev_prediction = current_prediction;
        prev_seq_id = ( *itr_seq)->GetSeqID();
      }
      if( going_up_str.length() || length_str.length())
      {
        BCL_MessageStd( "Topology: " + going_up_str);
        BCL_MessageStd( "Sizes   : " + length_str);
      }
    }

    //! @brief assign ss for the contiguous segment of residues given
    //! @param SEGMENT_AA_PHI_PSI aas, phi, psi angles in the segment
    //! @return vector of predicted SS types
    storage::Vector< biol::SSType> CIPhiPsi::CalculateSolubleSSForSegment
    (
      storage::Vector< t_SegmentStorage> &SEGMENT_AA_PHI_PSI,
      const double &MEM_THICKNESS
    ) const
    {
      // only consider right-handed alpha (type I) helices
      static const size_t s_min_helix_size( 5);
      const size_t segment_size( SEGMENT_AA_PHI_PSI.GetSize());
      linal::Vector< size_t> could_be_helix( segment_size, size_t( 0));
      linal::Vector< double> psih( segment_size + 8, util::GetUndefined< double>());
      linal::Vector< double> z_coord( segment_size, util::GetUndefined< double>());
      size_t helix_stretch_begin( util::GetUndefined< size_t>()), res_index( 0);
      for
      (
        storage::Vector< t_SegmentStorage>::const_iterator
          itr( SEGMENT_AA_PHI_PSI.Begin()),
          itr_end( SEGMENT_AA_PHI_PSI.End());
        itr != itr_end;
        ++itr, ++res_index
      )
      {
        const double phi( itr->Second().First());
        const double psi( itr->Second().Second());

        // determine whether this residue could be a helical residue
        if( phi < 180.0 || psi >= 45.0)
        {
          // could not be a helical residue
          // if this terminates a current helical segment, test its length
          if( util::IsDefined( helix_stretch_begin) && res_index - helix_stretch_begin >= s_min_helix_size)
          {
            // could be a helix, update could_be_helix
            for(; helix_stretch_begin < res_index; ++helix_stretch_begin)
            {
              could_be_helix( helix_stretch_begin) = size_t( 1);
            }
          }
          else
          {
            // non-helical; remove psi angles
            for( ; helix_stretch_begin < res_index; ++helix_stretch_begin)
            {
              psih( helix_stretch_begin + 4) = util::GetUndefined< double>();
            }
          }
          helix_stretch_begin = util::GetUndefined< size_t>();

          continue;
        }

        // could be a helical residue, update helix_stretch_begin if it is not currently defined
        if( !util::IsDefined( helix_stretch_begin))
        {
          helix_stretch_begin = res_index;
        }
        psih( res_index + 4) = psi;

        // get a z-coordinate for membrane proteins
        if( MEM_THICKNESS)
        {
          double z( itr->First()->GetAtom( biol::GetAtomTypes().N).GetCoordinates().Z());
          if( !util::IsDefined( z))
          {
            z = itr->First()->GetAtom( biol::GetAtomTypes().CA).GetCoordinates().Z();
            if( !util::IsDefined( z))
            {
              z = itr->First()->GetAtom( biol::GetAtomTypes().C).GetCoordinates().Z();
              if( !util::IsDefined( z))
              {
                z = itr->First()->GetCenter().Z();
              }
            }
          }
          z_coord( res_index) = z;
        }
      }
      // if this terminates a current helical segment, test its length
      if( util::IsDefined( helix_stretch_begin))
      {
        if( could_be_helix( helix_stretch_begin))
        {
          ++helix_stretch_begin;
        }
        if( res_index - helix_stretch_begin >= s_min_helix_size)
        {
          // could be a helix, update could_be_helix
          for(; helix_stretch_begin < res_index; ++helix_stretch_begin)
          {
            could_be_helix( helix_stretch_begin) = size_t( 1);
          }
        }
        else
        {
          // non-helical; remove psi angles
          for( ; helix_stretch_begin < res_index; ++helix_stretch_begin)
          {
            psih( helix_stretch_begin + 4) = util::GetUndefined< double>();
          }
        }
      }

      for( size_t i( 0); i < size_t( 4); ++i)
      {
        psih( i) = psih( 7 - i);
        psih( segment_size + 7 - i) = psih( segment_size + i);
      }

      // remove sloppy helices
      // handle helix-like regions with poor psi angles.  Specifically, look up to 4 residues away, take a (triangularly)
      // weighted average of the psi angles.  If the weighted average of psi angles is greater than -11.75 degrees,
      // the helix is too perturbed and should be removed
      for
      (
        size_t res_index( 0), res_index_offset( 4);
        res_index < segment_size;
        ++res_index_offset, ++res_index
      )
      {
        if( !could_be_helix( res_index))
        {
          continue;
        }
        math::RunningAverage< double> psi_ave;
        for( int offset( -3), max_offset( 4); offset < max_offset; ++offset)
        {
          const double local_psi( psih( res_index_offset + offset));
          if( util::IsDefined( local_psi))
          {
            psi_ave.AddWeightedObservation( local_psi, 4.0 - math::Absolute( offset));
          }
          else if( offset < 0)
          {
            psi_ave.Reset();
          }
          else
          {
            break;
          }
        }

        if( psi_ave.GetAverage() > -11.75)
        {
          // Membrane regions are often resolved at much lower resolution than the rest of the structure, while the membrane
          // environment strongly favors helices over coils, so a more soft cutoff is used for determining helices inside
          // membrane regions
          if( math::Absolute( z_coord( res_index)) > MEM_THICKNESS || psi_ave.GetAverage() > 20.0)
          {
            could_be_helix( res_index) = 0;
            psih( res_index_offset) = util::GetUndefined< double>();
          }
        }
      }
      res_index = 0;

      linal::Vector< size_t> could_be_strand( segment_size, size_t( 0));
      for
      (
        storage::Vector< t_SegmentStorage>::const_iterator
          itr( SEGMENT_AA_PHI_PSI.Begin()), itr_end( SEGMENT_AA_PHI_PSI.End());
        itr != itr_end;
        ++itr, ++res_index
      )
      {
        if( itr->Third() > 0.5)
        {
          could_be_strand( res_index) = size_t( 1);
        }
      }
      // handle the bizarre case that there are residues which could be strands as well as helices
      if( linal::ScalarProduct( could_be_helix, could_be_strand))
      {
        for( res_index = 0; res_index < segment_size; ++res_index)
        {
          if( could_be_helix( res_index) && could_be_strand( res_index))
          {
            could_be_helix( res_index) = could_be_strand( res_index) = 0;
          }
        }
      }

      // check for undersized helices and strands; helices must be at least 5 residues, strands must be at least 3
      for( res_index = 0; res_index < segment_size; ++res_index)
      {
        if( could_be_helix( res_index))
        {
          size_t helix_start( res_index);
          while( res_index < segment_size && could_be_helix( res_index))
          {
            ++res_index;
          }
          const size_t helix_size( res_index - helix_start);
          if( helix_size < s_min_helix_size)
          {
            std::fill( could_be_helix.Begin() + helix_start, could_be_helix.Begin() + res_index, size_t( 0));
          }
          --res_index;
        }
        else if( could_be_strand( res_index))
        {
          size_t strand_start( res_index);
          while( res_index < segment_size && could_be_strand( res_index))
          {
            ++res_index;
          }
          const size_t strand_size( res_index - strand_start);
          if( strand_size < size_t( 3))
          {
            std::fill( could_be_strand.Begin() + strand_start, could_be_strand.Begin() + res_index, size_t( 0));
          }
          --res_index;
        }
      }

      // assemble final ss types
      storage::Vector< biol::SSType> ss_types( segment_size, biol::GetSSTypes().COIL);
      for( res_index = 0; res_index < segment_size; ++res_index)
      {
        if( could_be_helix( res_index))
        {
          ss_types( res_index) = biol::GetSSTypes().HELIX;
        }
        else if( could_be_strand( res_index))
        {
          ss_types( res_index) = biol::GetSSTypes().STRAND;
        }
      }

      return ss_types;
    }

    //! @brief assign ss for the contiguous segment of residues given
    //! @param SEGMENT_AA_PHI_PSI aas, phi, psi angles in the segment
    //! @return vector of predicted SS types
    storage::Vector< CIPhiPsi::TMDirectionTypeEnum> CIPhiPsi::CalculateTMDirectionsForSegment
    (
      storage::Vector< t_SegmentStorage> &SEGMENT_AA_PHI_PSI,
      const double &MEM_THICKNESS
    ) const
    {
      const size_t segment_size( SEGMENT_AA_PHI_PSI.GetSize());
      storage::Vector< CIPhiPsi::TMDirectionTypeEnum> z_coord_ids( segment_size, TMDirectionTypeEnum( e_NonMembrane));
      if( !MEM_THICKNESS)
      {
        return z_coord_ids;
      }
      // only consider right-handed alpha (type I) helices
      linal::Vector< double> z_coord_a( segment_size, util::GetUndefined< double>());
      linal::Vector< double> z_coord_b( segment_size, util::GetUndefined< double>());
      linal::Vector< double> z_coord_mn( segment_size, util::GetUndefined< double>());
      linal::Vector< double> z_coord_mx( segment_size, util::GetUndefined< double>());
      linal::Vector< double> z_coord_mm( segment_size, util::GetUndefined< double>());
      linal::Vector< int> z_coord_trans( segment_size, 0);
      size_t res_index( 0);
      for
      (
        storage::Vector< t_SegmentStorage>::const_iterator
          itr( SEGMENT_AA_PHI_PSI.Begin()),
          itr_end( SEGMENT_AA_PHI_PSI.End());
        itr != itr_end;
        ++itr, ++res_index
      )
      {
        // get a z-coordinate for membrane proteins
        double z( itr->First()->GetAtom( biol::GetAtomTypes().N).GetCoordinates().Z());
        if( !util::IsDefined( z))
        {
          z = itr->First()->GetAtom( biol::GetAtomTypes().CA).GetCoordinates().Z();
          if( !util::IsDefined( z))
          {
            z = itr->First()->GetAtom( biol::GetAtomTypes().C).GetCoordinates().Z();
            if( !util::IsDefined( z))
            {
              z = itr->First()->GetCenter().Z();
            }
          }
        }
        z_coord_a( res_index) = z;
        z = itr->First()->GetAtom( biol::GetAtomTypes().C).GetCoordinates().Z();
        if( !util::IsDefined( z))
        {
          z = itr->First()->GetAtom( biol::GetAtomTypes().CA).GetCoordinates().Z();
          if( !util::IsDefined( z))
          {
            z = itr->First()->GetAtom( biol::GetAtomTypes().N).GetCoordinates().Z();
            if( !util::IsDefined( z))
            {
              z = itr->First()->GetCenter().Z();
            }
          }
        }
        z_coord_b( res_index) = z;
      }
      for( size_t i( 0); i < segment_size; ++i)
      {
        math::RunningMinMax< double> mm;
        for( size_t j( std::max( i, size_t( 2)) - 2), j_end( std::min( segment_size, i + 3)); j < j_end; ++j)
        {
          mm += z_coord_a( j);
          mm += z_coord_b( j);
        }
        z_coord_mn( i) = mm.GetMin();
        z_coord_mx( i) = mm.GetMax();
        z_coord_mm( i) = mm.GetMin() + mm.GetMax();
      }
      const size_t segment_size_mm( segment_size - 1);
      if( math::Absolute( 0.5 * ( z_coord_a( 0) + z_coord_b( 0))) < MEM_THICKNESS)
      {
        z_coord_trans( 0) = z_coord_mm( 1) > z_coord_mm( 0) ? 1 : ( z_coord_mm( 1) < z_coord_mm( 0) ? -1 : 0);
      }
      for( size_t i( 1); i < segment_size_mm; ++i)
      {
        if( math::Absolute( 0.5 * ( z_coord_a( i) + z_coord_b( i))) < MEM_THICKNESS)
        {
          if( z_coord_trans( i - 1) > 0)
          {
            z_coord_trans( i) = z_coord_mm( i) >= z_coord_mm( i - 1) && z_coord_mm( i + 1) >= z_coord_mm( i)
                                 ? 1
                                 : ( z_coord_mm( i) < z_coord_mm( i - 1) && z_coord_mm( i + 1) < z_coord_mm( i) ? -1 : 0);
            if( z_coord_trans( i) > 0)
            {
              z_coord_trans( i) += z_coord_trans( i - 1);
            }
          }
          else if( z_coord_trans( i - 1) < 0)
          {
            z_coord_trans( i) = z_coord_mm( i) > z_coord_mm( i - 1) && z_coord_mm( i + 1) > z_coord_mm( i)
                                 ? 1
                                 : ( z_coord_mm( i) <= z_coord_mm( i - 1) && z_coord_mm( i + 1) <= z_coord_mm( i) ? -1 : 0);
            if( z_coord_trans( i) && z_coord_trans( i) < 0)
            {
              z_coord_trans( i) += z_coord_trans( i - 1);
            }
          }
          else
          {
            z_coord_trans( i) = z_coord_mm( i) > z_coord_mm( i - 1) && z_coord_mm( i + 1) > z_coord_mm( i)
                                ? 1
                                : ( z_coord_mm( i) < z_coord_mm( i - 1) && z_coord_mm( i + 1) < z_coord_mm( i) ? -1 : 0);
          }
        }
      }
      if( math::Absolute( 0.5 * ( z_coord_a( segment_size_mm) + z_coord_b( segment_size_mm))) < MEM_THICKNESS)
      {
        z_coord_trans( segment_size_mm) = z_coord_trans( segment_size_mm - 1);
        z_coord_trans( segment_size_mm) += z_coord_trans( segment_size_mm) > 0
                                          ? 1
                                          : z_coord_trans( segment_size_mm) < 0 ? -1 : 0;
      }
      for( size_t i( segment_size - 1); i; --i)
      {
        if( z_coord_trans( i - 1) && z_coord_trans( i))
        {
          if( ( z_coord_trans( i) > 0) == ( z_coord_trans( i - 1) > 0))
          {
            z_coord_trans( i - 1) = z_coord_trans( i);
          }
        }
      }
      for( size_t i( 0); i < segment_size; ++i)
      {
        if( z_coord_trans( i) < 3 && z_coord_trans( i) > -3)
        {
          z_coord_trans( i) = 0;
        }
      }
      storage::Vector< storage::Triplet< TMDirectionTypeEnum, math::Range< size_t>, TMDirectionTypeEnum> > ranges;
      for( size_t i( 0); i < segment_size_mm; ++i)
      {
        if( z_coord_trans( i))
        {
          size_t i_start( i);
          math::RunningMinMax< double> z_mm;
          z_mm += z_coord_a( i);
          z_mm += z_coord_b( i);
          while( i < segment_size_mm && z_coord_trans( i) == z_coord_trans( i + 1))
          {
            i += 1;
            z_mm += z_coord_a( i);
            z_mm += z_coord_b( i);
          }
          TMDirectionTypeEnum label;
          if( z_mm.GetRange() < MEM_THICKNESS)
          {
            // insufficient vertical displacement. consider this to be a non-tm segment or an amphipathic region, depending
            // on length
            if( i - i_start + 1 > 6)
            {
              label = e_Amphipathic;
            }
            else
            {
              label = e_NonMembrane;
            }
          }
          else
          {
            label = z_coord_trans( i) > 0 ? e_LeavingCytosol : e_EnteringCytosol;
          }
          for( size_t j( i_start); j <= i; ++j)
          {
            z_coord_ids( j) = label;
          }
          if( label != e_NonMembrane)
          {
            ranges.PushBack
            (
              storage::Triplet< TMDirectionTypeEnum, math::Range< size_t>, TMDirectionTypeEnum>
              (
                label,
                math::Range< size_t>( i_start, i + 1),
                z_coord_trans( i) > 0 ? e_LeavingCytosol : e_EnteringCytosol
              )
            );
          }
        }
      }
      // fix rare cases where a helix that looks amphipathic is actually integral/effectively-transmembrane in terms of being
      // required to make the topology alternating as it should (e.g. Out-Amphi-Out should really be Out-In-Out)
      for( size_t i( 1), nr( std::max( ranges.GetSize(), size_t( 1)) - size_t( 1)); i < nr; ++i)
      {
        if
        (
          ranges( i).First() == e_Amphipathic
          && ranges( i - 1).First() != e_Amphipathic
          && ranges( i - 1).First() != ranges( i).Third()
          && ranges( i).Third() != ranges( i + 1).Third()
        )
        {
          ranges( i).First() = ranges( i).Third();
          TMDirectionTypeEnum label( ranges( i).Third());
          for( size_t j( ranges( i).Second().GetMin()), jm( ranges( i).Second().GetMax()); j < jm; ++j)
          {
            z_coord_ids( j) = label;
          }
        }
      }
      if( util::GetMessenger().GetCurrentMessageLevel() >= util::Message::e_Verbose)
      {
        std::string topo_str;
        size_t n_spans( 0);
        TMDirectionTypeEnum prev_direction( e_NonMembrane);
        size_t prev_direction_start( 0);
        std::string span_range_str( util::Format()( SEGMENT_AA_PHI_PSI( 0).First()->GetChainID()) + " : ");
        for( size_t i( 0); i < segment_size; ++i)
        {
          BCL_MessageStd
          (
            SEGMENT_AA_PHI_PSI( i).First()->GetIdentification()
            + " " + GetTMDirectionType( z_coord_ids( i)) + " "
            + " " + util::Format()( z_coord_mm( i)) + " "
            + util::Format()( z_coord_mn( i)) + " " + util::Format()( z_coord_mx( i))
            + " " + util::Format()( z_coord_trans( i))
          );
          topo_str += GetTMDirectionType( z_coord_ids( i));
          if( prev_direction != z_coord_ids( i))
          {
            if( prev_direction == e_EnteringCytosol || prev_direction == e_LeavingCytosol)
            {
              span_range_str += prev_direction.GetString()
                                + " " + util::Format()( SEGMENT_AA_PHI_PSI( prev_direction_start).First()->GetSeqID())
                                + " " + util::Format()( SEGMENT_AA_PHI_PSI( i - 1).First()->GetSeqID()) + " : ";
              ++n_spans;
            }
            prev_direction = z_coord_ids( i);
            prev_direction_start = i;
          }
        }
        if( prev_direction == e_EnteringCytosol || prev_direction == e_LeavingCytosol)
        {
          span_range_str += prev_direction.GetString()
                            + " " + util::Format()( SEGMENT_AA_PHI_PSI( prev_direction_start).First()->GetSeqID())
                            + " " + util::Format()( SEGMENT_AA_PHI_PSI( segment_size_mm).First()->GetSeqID()) + " : ";
          ++n_spans;
        }
        BCL_MessageStd( util::Format()( SEGMENT_AA_PHI_PSI( 0).First()->GetChainID()) + " topo_str " + topo_str);
        BCL_MessageStd( util::Format()( n_spans) + " spans: " + span_range_str);
      }
      return z_coord_ids;
    }

  } // namespace sspred
} // namespace bcl
