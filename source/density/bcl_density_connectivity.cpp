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
#include "density/bcl_density_connectivity.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically
#include <stack>

namespace bcl
{
  namespace density
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Connectivity::Connectivity()
    {
    }

    //! @brief Connectivity from various arguments
    Connectivity::Connectivity
    (
      const util::ShPtr< assemble::SSEGeometryInterface> &BODY_A,
      const bool BODY_A_BEGIN,
      const util::ShPtr< assemble::SSEGeometryInterface> &BODY_B,
      const bool BODY_B_BEGIN,
      const double CONNECTIVITY,
      const double DISTANCE
    ) :
      m_BodyA( BODY_A),
      m_BodyABegin( BODY_A_BEGIN),
      m_BodyB( BODY_B),
      m_BodyBBegin( BODY_B_BEGIN),
      m_Connectivity( CONNECTIVITY),
      m_Distance( DISTANCE)
    {
    }

    //! @brief determine all 4 Connectivities from every set of two bodies in a list in a density map
    //! @param DENSITY_MAP the density map the bodies are lying in
    //! @param BODIES list of bodies for which the mutual connectivities are to be determined
    //! @return list of connectivities of combined begin-begin, begin-end, end-begin, end-end for each pair of bodies
    storage::List< Connectivity>
    Connectivity::DetermineConnectivities
    (
      const Map &DENSITY_MAP,
      const util::ShPtrList< assemble::SSEGeometryInterface> &BODIES
    )
    {
      // list of all connectivities
      storage::List< Connectivity> all_connectivities;

      // iterate over all bodies
      for
      (
        util::ShPtrList< assemble::SSEGeometryInterface>::const_iterator body_itr1( BODIES.Begin()), body_itr_end( BODIES.End());
        body_itr1 != body_itr_end;
        ++body_itr1
      )
      {
        util::ShPtrList< assemble::SSEGeometryInterface>::const_iterator body_itr2( body_itr1);
        ++body_itr2;
        for
        (
          ;
          body_itr2 != body_itr_end;
          ++body_itr2
        )
        {
          // determine connectivities for that pair of bodies
          all_connectivities.Append( DetermineConnectivities( DENSITY_MAP, *body_itr1, *body_itr2));
        }
      }

      // end
      return all_connectivities;
    }

    //! @brief determine all 4 DensityConnectivities from 2 bodies in a density map given an optional length reduction
    //! @param DENSITY_MAP the density map the bodies are lying in
    //! @param BODY_A body a
    //! @param BODY_B body b
    //! @return list of connectivities of combined begin-begin, begin-end, end-begin, end-end of both bodies
    storage::List< Connectivity>
    Connectivity::DetermineConnectivities
    (
      const Map &DENSITY_MAP,
      const util::ShPtr< assemble::SSEGeometryInterface> &BODY_A,
      const util::ShPtr< assemble::SSEGeometryInterface> &BODY_B
    )
    {
      // list of all connectivities
      storage::List< Connectivity> all_connectivities;

      // find point at beginning of body_a to calculate connectivity (more shifted towards center of density rod)
      const linal::Vector3D body_a_begin_connectivity
      (
        // use center of first fragment (vs end point) because the first fragment is generally shifted into sse
        BODY_A->GetGeometries().FirstElement()->GetCenter()
//        BODY_A->GetGeometries().FirstElement()->GetMainAxis().GetEndPoint()
      );
      BCL_MessageDbg
      (
        "ATOM   9996  N   ASP B   1      " + util::Format().W( 6).FFP( 3).R()( body_a_begin_connectivity.X()) +
        "  " + util::Format().W( 6).FFP( 3).R()( body_a_begin_connectivity.Y()) +
        "  " + util::Format().W( 6).FFP( 3).R()( body_a_begin_connectivity.Z()) + "  1.00 38.08           N"
      );

      // find point at beginning of body_a to calculate distance (more shifted towards tip of density rod)
      const linal::Vector3D body_a_begin_distance
      (
        BODY_A->GetGeometries().FirstElement()->GetMainAxis().GetStartPoint()
      );

      // find point at end of body_a to calculate connectivity (more shifted towards center of density rod)
      const linal::Vector3D body_a_end_connectivity
      (
        BODY_A->GetGeometries().LastElement()->GetMainAxis().GetStartPoint()
      );
      BCL_MessageDbg
      (
        "ATOM   9996  N   ASP B   1      " + util::Format().W( 6).FFP( 3).R()( body_a_end_connectivity.X()) +
        "  " + util::Format().W( 6).FFP( 3).R()( body_a_end_connectivity.Y()) +
        "  " + util::Format().W( 6).FFP( 3).R()( body_a_end_connectivity.Z()) + "  1.00 38.08           N"
      );

      // find point at end of body_a to calculate distance (more shifted towards tip of density rod)
      const linal::Vector3D body_a_end_distance
      (
        BODY_A->GetGeometries().LastElement()->GetMainAxis().GetEndPoint()
      );

      // find point at beginning of body_b to calculate connectivity (more shifted towards center of density rod)
      const linal::Vector3D body_b_begin_connectivity
      (
        // use center of first fragment (vs end point) because the first fragment is generally shifted into sse
        BODY_B->GetGeometries().FirstElement()->GetCenter()
//        BODY_B->GetGeometries().FirstElement()->GetMainAxis().GetEndPoint()
      );
      BCL_MessageDbg
      (
        "ATOM   9996  N   ASP B   1      " + util::Format().W( 6).FFP( 3).R()( body_b_begin_connectivity.X()) +
        "  " + util::Format().W( 6).FFP( 3).R()( body_b_begin_connectivity.Y()) +
        "  " + util::Format().W( 6).FFP( 3).R()( body_b_begin_connectivity.Z()) + "  1.00 38.08           N"
      );

      // find point at beginning of body_b to calculate distance (more shifted towards tip of density rod)
      const linal::Vector3D body_b_begin_distance
      (
        BODY_B->GetGeometries().FirstElement()->GetMainAxis().GetStartPoint()
      );

      // find point at end of body_b to calculate connectivity (more shifted towards center of density rod)
      const linal::Vector3D body_b_end_connectivity
      (
        BODY_B->GetGeometries().LastElement()->GetMainAxis().GetStartPoint()
      );
      BCL_MessageDbg
      (
        "ATOM   9996  N   ASP B   1      " + util::Format().W( 6).FFP( 3).R()( body_b_end_connectivity.X()) +
        "  " + util::Format().W( 6).FFP( 3).R()( body_b_end_connectivity.Y()) +
        "  " + util::Format().W( 6).FFP( 3).R()( body_b_end_connectivity.Z()) + "  1.00 38.08           N"
      );

      // find point at end of body_b to calculate distance (more shifted towards tip of density rod)
      const linal::Vector3D body_b_end_distance
      (
        BODY_B->GetGeometries().LastElement()->GetMainAxis().GetEndPoint()
      );

      // create Connectivity object for each of the 4 possible connections between the bodies and insert them in
      // list of DensityConnectivities for return
      Connectivity begin_a_begin_b
      (
        BODY_A,
        true,
        BODY_B,
        true,
        ConnectionIntensity( DENSITY_MAP, body_a_begin_connectivity, body_b_begin_connectivity),
        linal::Distance( body_a_begin_distance, body_b_begin_distance)
      );
      all_connectivities.PushBack( std::move( begin_a_begin_b));

      Connectivity begin_a_end_b
      (
        BODY_A,
        true,
        BODY_B,
        false,
        ConnectionIntensity( DENSITY_MAP, body_a_begin_connectivity, body_b_end_connectivity),
        linal::Distance( body_a_begin_distance, body_b_end_distance)
      );
      all_connectivities.PushBack( std::move( begin_a_end_b));

      Connectivity end_a_begin_b
      (
        BODY_A,
        false,
        BODY_B,
        true,
        ConnectionIntensity( DENSITY_MAP, body_a_end_connectivity, body_b_begin_connectivity),
        linal::Distance( body_a_end_distance, body_b_begin_distance)
      );
      all_connectivities.PushBack( std::move( end_a_begin_b));

      Connectivity end_a_end_b
      (
        BODY_A,
        false,
        BODY_B,
        false,
        ConnectionIntensity( DENSITY_MAP, body_a_end_connectivity, body_b_end_connectivity),
        linal::Distance( body_a_end_distance, body_b_end_distance)
      );
      all_connectivities.PushBack( std::move( end_a_end_b));

      // end
      return all_connectivities;
    }

    //! @brief determine the lowest intensity the is necessary to connect two points in the electron density map
    //!        (this is the weakest point in the connection between two points)
    //! @param DENSITY_MAP the density map where the points are in
    //! @param COORDINATE_A the first point
    //! @param COORDINATE_B the second point
    //! @return an intensity that is low enough to connect the two points through voxels of intensities larger than the
    //!         returned intensity
    double
    Connectivity::ConnectionIntensity
    (
      const Map &DENSITY_MAP,
      const linal::Vector3D &COORDINATE_A,
      const linal::Vector3D &COORDINATE_B
    )
    {
      // retrieve a subdensity that contains the two points
      // (this is done to save time, reasoning is that the connection will appear in vicinity of both points)
      const Map sub_density( SubDensityFromTwoCoordinates( DENSITY_MAP, COORDINATE_A, COORDINATE_B));

      // retrieve the index of the two points within this subdensity
      linal::Vector< int> orig_index( 3);
      orig_index( 0) = sub_density.GetIndex()( 0);
      orig_index( 1) = sub_density.GetIndex()( 1);
      orig_index( 2) = sub_density.GetIndex()( 2);
      const linal::Vector< int> index_a_lnl( linal::Vector< int>( COORDINATE_A / sub_density.GetCellWidth()) - orig_index);
      const linal::Vector< int> index_b_lnl( linal::Vector< int>( COORDINATE_B / sub_density.GetCellWidth()) - orig_index);
      const linal::Vector< int> index_a( index_a_lnl.Begin(), index_a_lnl.End());
      const linal::Vector< int> index_b( index_b_lnl.Begin(), index_b_lnl.End());

      // check whether both points lie within subdensity map, if not return undefined intensity
      if
      (
        index_a( 0) < 0 || index_a( 0) >= int( sub_density.GetDimensions()( 0)) ||
        index_a( 1) < 0 || index_a( 1) >= int( sub_density.GetDimensions()( 1)) ||
        index_a( 2) < 0 || index_a( 2) >= int( sub_density.GetDimensions()( 2)) ||
        index_b( 0) < 0 || index_b( 0) >= int( sub_density.GetDimensions()( 0)) ||
        index_b( 1) < 0 || index_b( 1) >= int( sub_density.GetDimensions()( 1)) ||
        index_b( 2) < 0 || index_b( 2) >= int( sub_density.GetDimensions()( 2))
      )
      {
        BCL_MessageCrt( " indices outside of subdensity");
        return util::GetUndefined< double>();
      }

      // determine if a or b is the point with higher intensity
      linal::Vector< int> high_intensity_point( 3);
      linal::Vector< int> low_intensity_point( 3);

      if( sub_density( index_a( 0), index_a( 1), index_a( 2)) < sub_density( index_b( 0), index_b( 1), index_b( 2)))
      {
        high_intensity_point = index_b;
        low_intensity_point  = index_a;
      }
      else
      {
        high_intensity_point = index_a;
        low_intensity_point  = index_b;
      }

      // intensity at the lower intensity point
      const double low_intensity( sub_density( low_intensity_point( 0), low_intensity_point( 1), low_intensity_point( 2)));

      // stepsize by which intensity is lowered to check whether connection between the two points has already been made
      // this starts at low_intensity (have to include at least both points!) and goes down to the minimum intensity in
      // subdensity map
      const double stepsize( ( low_intensity - sub_density.GetMinimum()) / 100);

      // iterate over the intensities, this starts at low_intensity (have to include at least both points!)
      // and goes down to the minimum intensity in subdensity map
      for
      (
        double current_intensity = low_intensity; current_intensity >= sub_density.GetMinimum();
           current_intensity -= stepsize)
      {
        // make temp copy of density map (to color connected voxels as undefined)
        Map tmp_density( sub_density);

        // call iterative function ColorNeighboringVoxels, this will "color" all voxels as undefined that can be reached
        // from the high_intensity_point given the current intensity
        ColorNeighboringVoxels( tmp_density, current_intensity, high_intensity_point);

        // if current intensity was sufficiently low to connect both points the low_intensity_point will be colored as
        // undefined; thus return current_intensity in that case
        if( !util::IsDefined( tmp_density( low_intensity_point( 0), low_intensity_point( 1), low_intensity_point( 2))))
        {
          return current_intensity;
        }
      }

      // no connection was made
      return util::GetUndefined< double>();
    }

    //! @brief creates a subdensity that encapsulates the two coordinates in a 2 times the distance between the points large densitymap
    //! @param DENSITY_MAP the density map where the points are in
    //! @param COORDINATE_A the first point
    //! @param COORDINATE_B the second point
    //! @return the subdensity which contains the two points with generous spacing
    Map
    Connectivity::SubDensityFromTwoCoordinates
    (
      const Map &DENSITY_MAP,
      const linal::Vector3D &COORDINATE_A,
      const linal::Vector3D &COORDINATE_B
    )
    {
      // store the index for the lower right and upper left corner of the argument density map
      linal::Vector< int> lower_right_corner_index( 3, int( 0));
      linal::Vector< int> upper_left_corner_index( 3, int( 0));

      // determine the minimum and maximum index required
      for( size_t i( 0); i < 3; ++i)
      {
        lower_right_corner_index( i) = int( std::floor( std::min( COORDINATE_A( i), COORDINATE_B( i)) / DENSITY_MAP.GetCellWidth()( i))) - int( DENSITY_MAP.GetIndex()( i));
        upper_left_corner_index( i) = int( std::ceil( std::max( COORDINATE_A( i), COORDINATE_B( i)) / DENSITY_MAP.GetCellWidth()( i))) - int( DENSITY_MAP.GetIndex()( i));
      }

      // extension from the difference of the indices
      const linal::Vector< int> extension_subdensity( ( upper_left_corner_index - lower_right_corner_index) + 1);

      //return the sub density from the calculated parameters
      return DENSITY_MAP.SubMap
        (
          lower_right_corner_index( 0), lower_right_corner_index( 1), lower_right_corner_index( 2),
          extension_subdensity( 0), extension_subdensity( 1), extension_subdensity( 2)
        );
    }

    //! @brief set the voxel for the current point to "nan" if its intensity is larger than the threshold
    //! @param DENSITY_MAP the map that contains the voxel
    //! @param THRESHOLD threshold for the intensity that determines if the point should be colored or not
    //! @param CURRENT_POINT the index of the voxel to be checked
    //! @return False, if it was already colored or if it is smaller than the threshold. True if it was colored.
    bool
    Connectivity::ColorVoxel
    (
      Map &DENSITY_MAP,
      const double THRESHOLD,
      const linal::Vector< int> &CURRENT_POINT
    )
    {
      // calculate the intensity at the CURRENT_POINT
      double &current_intensity( DENSITY_MAP( CURRENT_POINT( 0), CURRENT_POINT( 1), CURRENT_POINT( 2)));

      // return false if it was already colored (or the intensity is undefined for some other reason)
      if( !util::IsDefined( current_intensity))
      {
        return false;
      }

      // if intensity of CURRENT_POINT is larger than the threshold return true and color the voxel (set it undefined)
      if( current_intensity >= THRESHOLD)
      {
        current_intensity = util::GetUndefined< double>();
        return true;
      }

      // return false if intensity of CURRENT_POINT is lower than threshold
      return false;
    }

    //! @brief iterate over all neighboring voxels of the MIDDLE_POINT
    //! @param DENSITY_MAP map where the MIDDLE_POINT lies in
    //! @param THRESHOLD threshold for the coloring (all voxels that have intensity THRESHOLD and higher are colored)
    //! @param MIDDLE_POINT the center point from which the surrounding voxels are explored and colored if above THRESHOLD
    inline
    void
    Connectivity::ColorNeighboringVoxels
    (
      Map &DENSITY_MAP,
      const double THRESHOLD,
      const linal::Vector< int> &MIDDLE_POINT
    )
    {
      // initialize empty stack that will hold all the voxels whose environment should be checked
      std::stack< linal::Vector< int> > voxels_to_be_checked;

      // add MIDDLE_POINT to stack
      voxels_to_be_checked.push( MIDDLE_POINT);

      // call ColorVoxel for the MIDDLE_POINT (this should always be colored)
      ColorVoxel( DENSITY_MAP, THRESHOLD, MIDDLE_POINT);

      // as long as there are voxels that need to be checked in the stack stay in while loop
      while( !voxels_to_be_checked.empty())
      {
        // remove top element from stack of voxels that need to checked and save it to voxel_to_be_checked
        linal::Vector< int> voxel_to_be_checked( voxels_to_be_checked.top());
        voxels_to_be_checked.pop();

        // initialize empty vector of surrounding voxels
        linal::Vector< int> environment_voxel_to_be_checked( 3);

        // iterate over the surrounding 27 voxels around voxel_to_be_checked
        // iterate over the surrounding voxels in the first dimension
        for( int i = voxel_to_be_checked( 0) - 1; i <= voxel_to_be_checked( 0) + 1; ++i)
        {
          // skip if you go outside the limits of the DENSITY_MAP
          if( i < 0 || i >= int( DENSITY_MAP.GetDimensions()( 0)))
          {
            continue;
          }
          // initialize first coordinate of environment_voxel_to_be_checked
          environment_voxel_to_be_checked( 0) = i;

          // iterate over the surrounding voxels in the second dimension
          for( int j = voxel_to_be_checked( 1) - 1; j <= voxel_to_be_checked( 1) + 1; ++j)
          {
            // skip if you go outside the limits of the DENSITY_MAP
            if( j < 0 || j >= int( DENSITY_MAP.GetDimensions()( 1)))
            {
              continue;
            }
            // initialize second coordinate of environment_voxel_to_be_checked
            environment_voxel_to_be_checked( 1) = j;

            // iterate over the surrounding voxels in the third dimension
            for( int k = voxel_to_be_checked( 2) - 1; k <= voxel_to_be_checked( 2) + 1; ++k)
            {
              // skip if you go outside the limits of the DENSITY_MAP
              if( k < 0 || k >= int( DENSITY_MAP.GetDimensions()( 2)))
              {
                continue;
              }
              // initialize third coordinate of environment_voxel_to_be_checked
              environment_voxel_to_be_checked( 2) = k;

              // color the environment voxel if it is below threshold and add it to stack if colored
              if( ColorVoxel( DENSITY_MAP, THRESHOLD, environment_voxel_to_be_checked))
              {
                voxels_to_be_checked.push( environment_voxel_to_be_checked);
              }
            }
          }
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read Connectivity from std::istream
    std::istream &Connectivity::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_BodyA, ISTREAM);
      io::Serialize::Read( m_BodyABegin, ISTREAM);
      io::Serialize::Read( m_BodyB, ISTREAM);
      io::Serialize::Read( m_BodyBBegin, ISTREAM);
      io::Serialize::Read( m_Connectivity, ISTREAM);
      io::Serialize::Read( m_Distance, ISTREAM);

      // end
      return ISTREAM;
    }

    //! write Connectivity into std::ostream
    std::ostream &Connectivity::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // read member
      io::Serialize::Write( m_BodyA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BodyABegin, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BodyB, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BodyBBegin, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Connectivity, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Distance, OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

    //! @brief boolean operator CONNECTIVITY_LHS == CONNECTIVITY_RHS
    //! @param CONNECTIVITY_LHS first CONNECTIVITY
    //! @param CONNECTIVITY_RHS second CONNECTIVITY
    //! @return whether CONNECTIVITY_LHS is equal to CONNECTIVITY_RHS
    bool operator ==( const Connectivity &CONNECTIVITY_LHS, const Connectivity &CONNECTIVITY_RHS)
    {
      return
      (
        CONNECTIVITY_LHS.GetBodyA() == CONNECTIVITY_RHS.GetBodyA() &&
        CONNECTIVITY_LHS.GetBodyB() == CONNECTIVITY_RHS.GetBodyB() &&
        CONNECTIVITY_LHS.GetOrientationA() == CONNECTIVITY_RHS.GetOrientationA() &&
        CONNECTIVITY_LHS.GetOrientationB() == CONNECTIVITY_RHS.GetOrientationB() &&
        CONNECTIVITY_LHS.GetConnectivity() == CONNECTIVITY_RHS.GetConnectivity() &&
        CONNECTIVITY_LHS.GetDistance() == CONNECTIVITY_RHS.GetDistance()
      );
    }

  } // namespace density
} // namespace bcl
