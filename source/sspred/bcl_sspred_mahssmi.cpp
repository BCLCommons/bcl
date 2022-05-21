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
#include "sspred/bcl_sspred_mahssmi.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_classes.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "coord/bcl_coord_polygon.h"
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_const_graph.h"
#include "linal/bcl_linal_vector_2d_operations.h"
#include "math/bcl_math_limits.h"
#include "math/bcl_math_running_average.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "sspred/bcl_sspred_pdb.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Mahssmi::Mahssmi() :
      m_SSType(),
      m_EnvironmentType(),
      m_BetaBarrelResidueFacesPore( false),
      m_OriginatesInCytosol( false)
    {
    }

    //! @brief constructor from members
    //! @param SS_TYPE type of sse
    //! @param ENVIRONMENT environment of the residue
    //! @param FACES_PORE true if the residue faces the pore of the membrane
    //! @param ORIGINATES_CYTOSOL true if the residue is on a TM region that originates in the cytosol.
    //!        False for soluble regions and TM-regions that originate outside the cytosol
    Mahssmi::Mahssmi( const biol::SSType &SS_TYPE, const biol::EnvironmentType &ENVIRONMENT, const bool &FACES_PORE, const bool &ORIGINATES_CYTOSOL) :
      m_SSType( SS_TYPE),
      m_EnvironmentType( ENVIRONMENT),
      m_BetaBarrelResidueFacesPore( FACES_PORE),
      m_OriginatesInCytosol( ORIGINATES_CYTOSOL)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Mahssmi
    Mahssmi *Mahssmi::Clone() const
    {
      return new Mahssmi( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Mahssmi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &Mahssmi::GetFileExtension() const
    {
      static const std::string s_file_extension( ".mahssmi");

      // end
      return s_file_extension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get whether the residue faces the pore
    //! @return true if the residue is in a beta barrel and faces the pore
    bool Mahssmi::DoesBetaBarrelResidueFacePore() const
    {
      return m_BetaBarrelResidueFacesPore;
    }

    //! @brief get whether the residue is in a TM-region whose N-term is towards the cytosol
    //! @return true if the residue is in a TM-region whose N-term is towards the cytosol
    bool Mahssmi::OriginatesInCytosol() const
    {
      return m_OriginatesInCytosol;
    }

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D Mahssmi::GetThreeStatePrediction() const
    {
      return m_SSType->GetThreeStatePrediction();
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> Mahssmi::GetNineStatePrediction() const
    {
      return ConvertThreeStateToNineState( GetThreeStatePrediction(), m_EnvironmentType);
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &Mahssmi::ReadPredictionsForAA( std::istream &ISTREAM, biol::AABase &AMINO_ACID) const
    {
      std::string temp;

      // read chain id
      char chain_id;
      ISTREAM >> chain_id;
      while( ISTREAM.good() && chain_id != AMINO_ACID.GetChainID() && ( chain_id != '_' || AMINO_ACID.GetChainID() != ' '))
      {
        // call getline to read the next line
        std::getline( ISTREAM, temp);
        ISTREAM >> chain_id;
      }
      BCL_Assert
      (
        chain_id == AMINO_ACID.GetChainID() || ( AMINO_ACID.GetChainID() == ' ' && chain_id == '_'),
        "Chain ID mismatch; MAHSSMI file had " + util::Format()( chain_id)
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
        BCL_MessageStd( "Warning: Old MAHSMMI file used. This is only valid if this is not a membrane protein!");
      }
      bool originates_cytosol( temp[ 0] == 'c');

      AMINO_ACID.SetSSPrediction( GetMethods().e_MAHSSMI, Mahssmi( ss, env_type, faces_pore, originates_cytosol));

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &Mahssmi::ReadPredictionsForAASequence( std::istream &ISTREAM, biol::AASequence &AA_SEQUENCE) const
    {
      // remove header lines
      util::ChopHeader( ISTREAM);

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
        actual_chain_id == desired_chain_id,
        "Requested chain " + util::Format()( desired_chain_id) + " was not in MAHSSMI file!"
      );

      // get current aa type from read one letter code
      const biol::AAType current_aa_type( biol::GetAATypes().AATypeFromOneLetterCode( one_letter_code[ 0]));

      // get the first aa out of the sequence
      biol::AABase &first_aa( **AA_SEQUENCE.Begin());

      // assert matching amino acid seqids
      BCL_Assert
      (
        first_aa.GetSeqID() == seqid,
        "mismatch in seqids!\n sequence: " + first_aa.GetIdentification() +
        " vs. from sspred MAHSSMI: " + util::Format()( seqid) + " " + one_letter_code
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
        MethodHandler::ReadPredictionsForAA( ISTREAM, **aa_itr, GetMethods().e_MAHSSMI);
      }

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

    //! @brief iterates over the sequences in ProteinModel and calculates MAHSSMI for every residue in the sequence
    //! @param PROTEIN_MODEL ProteinModel for which MAHSSMI will be calculated
    void Mahssmi::Calculate( assemble::ProteinModel &PROTEIN_MODEL, const bool &USE_PDBTM_MEMBRANE_THICKNESS)
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

    //! @brief iterates over the sequence and calculates MAHSSMI for every residue in the sequence
    //! @param SEQUENCE Sequence of interest
    //! @param MEMBRANE_PTR pointer to membrane object
    void Mahssmi::Calculate( biol::AASequence &SEQUENCE, const util::SiPtr< const biol::Membrane> &MEMBRANE_PTR)
    {
      // iterate over all residues in the given sequence
      for
      (
        biol::AASequence::iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // test whether this residue has a stride dssp identification
        util::SiPtr< const MethodInterface> dsspstride_id( ( *aa_itr)->GetSSPrediction( GetMethods().e_StrideDSSP));
        BCL_Assert
        (
          dsspstride_id.IsDefined(),
          "Missing DSSPStride prediction for residue " + ( *aa_itr)->GetIdentification()
        );

        // stride dssp prediction will be used for all soluble residues
        ( *aa_itr)->SetSSPrediction
        (
          GetMethods().e_MAHSSMI,
          Mahssmi
          (
            dsspstride_id->GetOneStateSSPrediction(),
            biol::GetEnvironmentTypes().e_Solution,
            false
          )
        );

        BCL_Assert
        (
          ( *aa_itr)->GetAAClass() == biol::GetAAClasses().e_AACaCb
          || ( *aa_itr)->GetAAClass() == biol::GetAAClasses().e_AAComplete
          || ( *aa_itr)->GetAAClass() == biol::GetAAClasses().e_AABackBone,
          "Computation of MAHSSMI requires -aaclass AACaCb or -aaclass AAComplete or -aaclass AABackBone, but was "
          + ( *aa_itr)->GetAAClass().GetName()
        );
      }

      if( !MEMBRANE_PTR.IsDefined())
      {
        return;
      }

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
      size_t tm_strand_count( 0), tm_helix_count( 0);

      // get the palsse types
      // use the palsse identified values only in the membrane

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

        // test whether this residue has a palsse identification
        util::SiPtr< const MethodInterface> palsse_id( ( *aa_itr)->GetSSPrediction( GetMethods().e_PALSSE));
        BCL_Assert
        (
          palsse_id.IsDefined(),
          "Missing palsse id for residue " + ( *aa_itr)->GetIdentification()
        );

        // mahssmi only differs from stridedssp in the membrane core region, so skip residues that do not appear to be
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
        const biol::SSType ss( palsse_id->GetOneStateSSPrediction());

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

        // update the MAHSSMI analysis to denote that this is a membrane core residue
        // This will only be changed if the residue is detected to be in a pore
        ( *aa_itr)->SetSSPrediction
        (
          GetMethods().e_MAHSSMI,
          Mahssmi( ss, biol::GetEnvironmentTypes().e_MembraneCore, false)
        );

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
        ComputeCytosolicTMOrientation( SEQUENCE);
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
          GetMethods().e_MAHSSMI,
          Mahssmi
          (
            ( *aa_itr)->GetSSPrediction( GetMethods().e_StrideDSSP)->GetOneStateSSPrediction(),
            biol::GetEnvironmentTypes().e_Solution,
            false
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
              GetMethods().e_MAHSSMI,
              Mahssmi
              (
                ( *itr_aa)->GetSSPrediction( GetMethods().e_StrideDSSP)->GetOneStateSSPrediction(),
                biol::GetEnvironmentTypes().e_Solution,
                false
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
          if( is_closer)
          {
            ( *itr_aa)->SetSSPrediction
            (
              GetMethods().e_MAHSSMI,
              Mahssmi
              (
                ( *itr_aa)->GetSSPrediction( GetMethods().e_PALSSE)->GetOneStateSSPrediction(),
                biol::GetEnvironmentTypes().e_MembraneCore,
                true
              )
            );
          }
        }
      }

      // for each residue in a membrane strand
      // iterate over all residues in the given sequence
      for
      (
        biol::AASequence::iterator
          aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // get a const-ref to the coordinates of this residues' CA
        const linal::Vector3D &coordinates( ( *aa_itr)->GetCA().GetCoordinates());

        // determine the z-coordinate of the previous aas
        const double z_coordinate( coordinates.Z());
        const double z_coordinate_abs( math::Absolute( z_coordinate));
        if( z_coordinate_abs > membrane_thickness)
        {
          // skip residues that have to be outside the membrane
          continue;
        }

        util::SiPtr< const MethodInterface> mahssmi_si( ( *aa_itr)->GetSSPrediction( GetMethods().e_MAHSSMI));
        // skip non-membrane residues
        if
        (
          !mahssmi_si.IsDefined()
          || mahssmi_si->GetOneStateTMPrediction() != biol::GetEnvironmentTypes().e_MembraneCore
        )
        {
          continue;
        }

        const biol::SSType ss_type( mahssmi_si->GetOneStateSSPrediction());

        // extract up to the next 4 residues ss types, z-coordinates, and is-in-membrane
        biol::AASequence::iterator aa_itr_next( aa_itr);
        ++aa_itr_next;
        std::string next_sstypes;
        std::string next_membrane;
        storage::Vector< double> next_z_coordinates;
        storage::Vector< double> next_z_coordinates_abs;
        next_sstypes += ss_type->GetOneLetterCode();
        next_membrane += 'M';
        next_z_coordinates.PushBack( z_coordinate);
        next_z_coordinates_abs.PushBack( z_coordinate_abs);
        const int seqid( ( *aa_itr)->GetSeqID());
        bool missing_coordinates( false);
        for
        (
          size_t n_residues_extracted( 0);
          aa_itr_next != aa_itr_end;
          ++n_residues_extracted, ++aa_itr_next
        )
        {
          // ensure that seq ids are continuous
          if( ( *aa_itr_next)->GetSeqID() != seqid + 1 + int( n_residues_extracted))
          {
            missing_coordinates = true;
            break;
          }

          util::SiPtr< const MethodInterface> mahssmi_si_next
          (
            ( *aa_itr_next)->GetSSPrediction( GetMethods().e_MAHSSMI)
          );
          if( !mahssmi_si_next.IsDefined())
          {
            missing_coordinates = true;
            break;
          }
          // get a const-ref to the coordinates of this residues' CA
          const double &next_coordinate_z( ( *aa_itr_next)->GetCA().GetCoordinates().Z());
          if( !util::IsDefined( next_coordinate_z))
          {
            missing_coordinates = true;
            break;
          }
          // stop on reaching a helix
          const char onelettercode( mahssmi_si_next->GetOneStateSSPrediction()->GetOneLetterCode());
          if( onelettercode == 'H')
          {
            break;
          }
          if( math::Absolute( next_coordinate_z) >= membrane_thickness)
          {
            break;
          }

          // stop if the next residue is clearly makes a turn
          if( n_residues_extracted > 1)
          {
            const bool dominant_direction( next_z_coordinates( n_residues_extracted - ( n_residues_extracted & 1)) > next_z_coordinates( 0));
            if( ( next_coordinate_z > next_z_coordinates( n_residues_extracted - 1)) != dominant_direction)
            {
              break;
            }
          }
          next_membrane +=
            mahssmi_si_next->GetOneStateTMPrediction() == biol::GetEnvironmentTypes().e_MembraneCore ? 'M' : 'S';
          next_sstypes += onelettercode;

          // determine the z-coordinate of the previous aas
          next_z_coordinates.PushBack( next_coordinate_z);
          next_z_coordinates_abs.PushBack( math::Absolute( next_coordinate_z));
        }

        // if we're missing coordinates, there are no heuristics based on slope and size that can rationally be used
        if( missing_coordinates)
        {
          continue;
        }

        // trim soluble residues at the end
        size_t n_residues_in_stretch( next_membrane.rfind( 'M') + 1);
        if( n_residues_in_stretch != next_membrane.size())
        {
          const size_t n_to_erase( next_membrane.size() - n_residues_in_stretch);
          next_membrane.erase( n_residues_in_stretch, n_to_erase);
          next_sstypes.erase( n_residues_in_stretch, n_to_erase);
          next_z_coordinates.RemoveElements( n_residues_in_stretch, n_to_erase);
          next_z_coordinates_abs.RemoveElements( n_residues_in_stretch, n_to_erase);
        }

        linal::Vector< double> next_z_coords( next_z_coordinates);
        const double mn( next_z_coords.Min()), mx( next_z_coords.Max());
        const double slope( ( mx - mn) / ( std::max( n_residues_in_stretch, size_t( 2)) / 2));
        BCL_MessageVrb
        (
          "N-res: " + util::Format()( n_residues_in_stretch)
          + " min: " + util::Format()( mn)
          + " max: " + util::Format()( mx)
          + " span: " + util::Format()( mx - mn)
          + " slope: " + util::Format()( slope)
          + " start aa: " + ( *aa_itr)->GetIdentification()
          + " SS-types: " + next_sstypes
          + " Membrane: " + next_membrane
        );

        if
        (
          n_residues_in_stretch <= size_t( 2)
          ||
          ( n_residues_in_stretch == size_t( 3) && math::Absolute( next_z_coordinates( 0) - next_z_coordinates( 2)) < 2.5)
          ||
          ( slope < 2.5 && n_residues_in_stretch < 5)
        )
        {
          for( size_t n_residues_changed( 0); n_residues_changed < n_residues_in_stretch; ++aa_itr, ++n_residues_changed)
          {
            const biol::SSType ssprediction
            (
              ( *aa_itr)->GetSSPrediction( GetMethods().e_StrideDSSP)->GetOneStateSSPrediction()
            );
            // 1-2 residue in the membrane stretch; likely spurious, so remove them
            BCL_MessageVrb( ( *aa_itr)->GetIdentification() + " too short " + util::Format()( n_residues_in_stretch));
            ( *aa_itr)->SetSSPrediction
            (
              GetMethods().e_MAHSSMI,
              Mahssmi
              (
                ssprediction == biol::GetSSTypes().STRAND ? biol::GetSSTypes().COIL : ssprediction,
                biol::GetEnvironmentTypes().e_MembraneCore,
                false
              )
            );
          }
          --aa_itr;
          continue;
        }
        if( next_membrane[ 1] == 'S')
        {
          const size_t first_good_position( std::min( next_membrane.find( "MM"), n_residues_in_stretch));
          for( size_t pos( 0); pos < first_good_position; ++pos, ++aa_itr)
          {
            if( next_membrane[ pos] == 'M')
            {
              ( *aa_itr)->SetSSPrediction
              (
                GetMethods().e_MAHSSMI,
                Mahssmi
                (
                  ( *aa_itr)->GetSSPrediction( GetMethods().e_StrideDSSP)->GetOneStateSSPrediction(),
                  biol::GetEnvironmentTypes().e_Solution,
                  false
                )
              );
            }
          }
          --aa_itr;
          continue;
        }
        if( next_membrane[ n_residues_in_stretch - 2] == 'S')
        {
          for
          (
            size_t first_bad_position( next_membrane.rfind( "MM") + 3);
            first_bad_position < n_residues_in_stretch;
            ++first_bad_position
          )
          {
            if( next_membrane[ first_bad_position] == 'M')
            {
              biol::AASequence::iterator itr_aa( aa_itr + first_bad_position);
              ( *itr_aa)->SetSSPrediction
              (
                GetMethods().e_MAHSSMI,
                Mahssmi
                (
                  ( *itr_aa)->GetSSPrediction( GetMethods().e_StrideDSSP)->GetOneStateSSPrediction(),
                  biol::GetEnvironmentTypes().e_Solution,
                  false
                )
              );
            }
          }
          --aa_itr;
          continue;
        }

        // find the first strand residue in the set
        const size_t first_strand_res( next_sstypes.find( 'E')), last_strand_res( next_sstypes.rfind( 'E'));

        for( size_t pos( 0); pos < n_residues_in_stretch; ++pos, ++aa_itr)
        {
          const bool in_inner_strand( pos >= first_strand_res && pos <= last_strand_res);
          if( next_membrane[ pos] == 'S' || ( next_sstypes[ pos] == 'C' && in_inner_strand))
          {
            ( *aa_itr)->SetSSPrediction
            (
              GetMethods().e_MAHSSMI,
              Mahssmi
              (
                in_inner_strand ? biol::GetSSTypes().STRAND : biol::GetSSTypes().COIL,
                biol::GetEnvironmentTypes().e_MembraneCore,
                true
              )
            );
          }
        }
        if
        (
          aa_itr != aa_itr_end
          && ( *aa_itr)->GetSSPrediction( GetMethods().e_MAHSSMI)->GetOneStateSSPrediction() == biol::GetSSTypes().STRAND
          && ( *aa_itr)->GetSSPrediction( GetMethods().e_MAHSSMI)->GetOneStateTMPrediction() == biol::GetEnvironmentTypes().e_MembraneCore
          && next_sstypes[ next_sstypes.size() - 1] == 'E'
        )
        {
          --aa_itr;
          ( *aa_itr)->SetSSPrediction
          (
            GetMethods().e_MAHSSMI,
            Mahssmi
            (
              biol::GetSSTypes().COIL,
              biol::GetEnvironmentTypes().e_MembraneCore,
              true
            )
          );
        }
        else
        {
          --aa_itr;
        }
      }
      ComputeCytosolicTMOrientation( SEQUENCE);
    }

    //! @brief iterates over the sequence and writes MAHSSMI for every residue in every chain
    //! @param STREAM stream to write the analysis to
    //! @param PROTEIN_MODEL Model to write the analysis for
    void Mahssmi::WriteAnalysis( std::ostream &STREAM, const assemble::ProteinModel &PROTEIN_MODEL)
    {
      // write the normal header
      STREAM <<
       "# Membrane aware hybrid secondary structure and membrane identification file\n"
       "# Version 1.0\n"
       "# SecondaryStructure is one of H/E/C ( Helix / Strand / Coil )\n"
       "# MembraneLocalization is one of M/S ( Membrane / Solution )\n"
       "# Orientation is either L (residue faces lipid), P (residue faces pore), or -\n"
       "# Orientation of '-' is used for soluble residues\n"
       "# SeqID 1LtrCode Chain SecondaryStructure MembraneLocalization Orientation\n";

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

    //! @brief iterates over the protein and calculates MAHSSMI for every residue in every chain
    //! @param STREAM stream to write the analysis to
    //! @param SEQUENCE Sequence to write the analysis for
    //! @param CHAIN_ID Id of the chain that the sequence represents
    void Mahssmi::WriteAnalysis( std::ostream &STREAM, const biol::AASequence &SEQUENCE, const char &CHAIN_ID)
    {
      const char chain_id( CHAIN_ID == ' ' ? '_' : CHAIN_ID);
      // iterate over amino acids in the sequence
      for
      (
        biol::AASequence::const_iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        util::SiPtr< const Mahssmi> mahssmi_ptr( ( *aa_itr)->GetSSPrediction( GetMethods().e_MAHSSMI));

        // check existence of requested method
        if( !mahssmi_ptr.IsDefined())
        {
          BCL_MessageCrt( "MAHSSMI was not stored for " + ( *aa_itr)->GetIdentification());
          continue;
        }

        const Mahssmi &mahssmi( *mahssmi_ptr);

        // write residue info
        STREAM << ( *aa_itr)->GetSeqID() << ' ' << ( *aa_itr)->GetType()->GetOneLetterCode() << ' ' << chain_id << ' ';

        // write analysis info
        STREAM << mahssmi.GetOneStateSSPrediction()->GetOneLetterCode() << ' ';
        if( mahssmi.GetOneStateTMPrediction() == biol::GetEnvironmentTypes().e_MembraneCore)
        {
          STREAM << "M ";
          STREAM << ( mahssmi.DoesBetaBarrelResidueFacePore() ? 'P' : 'L') << ' ';
          STREAM << ( mahssmi.GetOneStateSSPrediction()->IsStructured() ? ( mahssmi.m_OriginatesInCytosol ? 'c' : 'e') : '-');
        }
        else
        {
          // soluble residue
          STREAM << "S - -";
        }
        STREAM << '\n';
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Mahssmi::Read( std::istream &ISTREAM)
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
    std::ostream &Mahssmi::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_EnvironmentType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BetaBarrelResidueFacesPore, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

    //! @brief set SSE Numbers for a given AA Sequence
    void Mahssmi::ComputeCytosolicTMOrientation( biol::AASequence &AA_SEQUENCE)
    {
      const Method mahssmi( GetMethods().e_MAHSSMI);
      biol::SSType prev_prediction( biol::GetSSTypes().COIL);
      biol::SSType current_prediction( biol::GetSSTypes().COIL);
      bool was_in_membrane( false);
      bool going_up( false);
      for
      (
        biol::AASequence::iterator itr_seq( AA_SEQUENCE.Begin()), itr_seq_end( AA_SEQUENCE.End());
        itr_seq != itr_seq_end;
        ++itr_seq
      )
      {
        util::SiPtr< const Mahssmi> sstype( ( *itr_seq)->GetSSPrediction( mahssmi));
        if( !sstype.IsDefined())
        {
          continue;
        }
        current_prediction = sstype->GetOneStateSSPrediction();
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
        if( in_membrane && !was_in_membrane)
        {
          going_up = false;
          // look ahead to the end of the membrane segment. If it is above this point, set going up flag
          biol::AASequence::iterator itr_seq_b( itr_seq);
          for( ; itr_seq_b != itr_seq_end; ++itr_seq_b)
          {
            util::SiPtr< const Mahssmi> sstype_b( ( *itr_seq_b)->GetSSPrediction( mahssmi));
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
          if
          (
            ( itr_seq_b - itr_seq) >= int( 3)
            && ( *itr_seq_b)->GetCA().GetCoordinates().Z() > ( *itr_seq)->GetCA().GetCoordinates().Z()
          )
          {
            going_up = true;
          }
        }
        was_in_membrane = in_membrane;
        Mahssmi newsstype( *sstype);
        newsstype.m_OriginatesInCytosol = going_up && in_membrane;
        ( *itr_seq)->SetSSPrediction( mahssmi, newsstype);

        util::SiPtr< const Mahssmi> sstypen( ( *itr_seq)->GetSSPrediction( mahssmi));
        prev_prediction = current_prediction;
      }
    }

  } // namespace sspred
} // namespace bcl
