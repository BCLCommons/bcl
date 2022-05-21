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
#include "density/bcl_density.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace density
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
#include "density/bcl_density_fit_protein_minimizer_mc.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "coord/bcl_coord_move_rotate_random.h"
#include "coord/bcl_coord_move_transform_random.h"
#include "coord/bcl_coord_move_translate_random.h"
#include "density/bcl_density_map.h"
#include "fold/bcl_fold_mutate_protein_model.h"
#include "math/bcl_math_mutate_decision_node.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_unimproved.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FitProteinMinimizerMC::s_Instance
    (
      GetObjectInstances().AddInstance( new FitProteinMinimizerMC())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FitProteinMinimizerMC::FitProteinMinimizerMC()
    {
    }

    //! @brief Clone function
    //! @return pointer to new FitProteinMinimizerMC
    FitProteinMinimizerMC *FitProteinMinimizerMC::Clone() const
    {
      return new FitProteinMinimizerMC( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FitProteinMinimizerMC::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the resolution of the density map
    //! @param RESOLUTION density map and simulation resolution
    void FitProteinMinimizerMC::SetResolution( const double RESOLUTION)
    {
      m_Resolution = RESOLUTION;
    }

    //! @brief set max translation and rotation
    //! @param MAX_TRANSLATION max translation in any direction for a single iteration
    //! @param MAX_ROTATION max rotation in radians in any direction for a single iteration
    void FitProteinMinimizerMC::SetMaxTranslationAndRotation( const double MAX_TRANSLATION, const double MAX_ROTATION)
    {
      m_MaxTranslation = MAX_TRANSLATION;
      m_MaxRotation    = MAX_ROTATION;
    }

    //! @brief set the max number of iterations for minimization
    //! @param MAX_NUMBER_ITERATIONS maximum number of iterations for minimization
    void FitProteinMinimizerMC::SetMaxIterations( const size_t MAX_NUMBER_ITERATIONS)
    {
      m_MaxIterations = MAX_NUMBER_ITERATIONS;
    }

    //! @brief set protein agreement measure to be used
    //! @param AGREEMENT protein agreement enumerator
    void FitProteinMinimizerMC::SetProteinAgreement( const ProteinAgreement &AGREEMENT)
    {
      m_ProteinAgreement = AGREEMENT;
    }

    //! @brief simulator to use
    //! @param DENSITY_SIMULATOR simulator enumerator
    void FitProteinMinimizerMC::SetSimulator( const Simulator &DENSITY_SIMULATOR)
    {
      m_Simulator = DENSITY_SIMULATOR;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator minimizing the position of a protein model within a given density map
    //! @param PROTEIN_MODEL start position of given protein model
    //! @param DENSITY_MAP the density map to fit the PROTEIN_MODEL into
    //! @return the fitted protein model
    assemble::ProteinModel FitProteinMinimizerMC::operator()( const assemble::ProteinModel &PROTEIN_MODEL, const Map &DENSITY_MAP) const
    {
      // construct a function that calculates the deviation between simulated density from atoms and the given density map
      util::ShPtr< ProteinAgreementInterface> density_agreement
      (
        GetProteinAgreements().CreateProteinAgreement
        (
          m_ProteinAgreement,
          m_Simulator,
          util::ToSiPtr( DENSITY_MAP),
          m_Resolution
        )
      );

      // mutate rotate
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > sp_mutate_rotate
      (
        new fold::MutateProteinModel( coord::MoveRotateRandom( m_MaxRotation))
      );

      // mutate translate
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > sp_mutate_translate
      (
        new fold::MutateProteinModel( coord::MoveTranslateRandom( m_MaxTranslation))
      );

      // mutate function rot trans
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > sp_mutate_rot_trans
      (
        new fold::MutateProteinModel( coord::MoveTransformRandom( m_MaxTranslation, m_MaxRotation))
      );

      math::ObjectProbabilityDistribution< math::MutateInterface< assemble::ProteinModel> > mutates;
      mutates.PushBack( 0.25, *sp_mutate_translate);
      mutates.PushBack( 0.25, *sp_mutate_rotate);
      mutates.PushBack( 0.25, *sp_mutate_rot_trans);

      // mutate function
      math::MutateDecisionNode< assemble::ProteinModel> mutate( mutates);

//      // create printer
//      util::ShPtr< mc::MoviePrinterInterface> sp_movie
//      (
//        new mc::MoviePrinterChimera()
//      );
//
//      // if user wished to print all intermediate minimization steps
//      if( m_WriteMinimizationFlag->GetFlag())
//      {
//        const storage::Set< opti::StepStatusEnum> step_status_set
//        (
//          opti::e_Improved, opti::e_Accepted, opti::e_Skipped, opti::e_Rejected
//        );
//
//        // initialize movie
//        sp_movie->Initialize
//        (
//          m_OutputPrefix->GetFirstParameter()->GetValue(),
//          storage::Vector< std::string>::Create( GetStaticClassName< storage::Table< double> >(), "value"),
//          storage::Vector< std::string>( 1, density_agreement->GetScheme()),
//          720, 480, false
//        );
//
//        // construct movie printer
//        util::ShPtr< assemble::PrinterProteinModelMovie> sp_printer
//        (
//          new assemble::PrinterProteinModelMovie
//          (
//            m_OutputPrefix->GetFirstParameter()->GetValue(),
//            sp_movie,
//            density_agreement,
//            step_status_set,
//            quality::GetSuperimposeMeasures().e_NoSuperimpose,
//            storage::Set< quality::Measure>()
//          )
//        );
//
//        // insert density map as first frame
//        sp_movie->SetStartFrame( m_MrcFilenameParam->GetValue(), true);
//
//        // set the water mark
//        sp_movie->SetWaterMark();
//
//        // set the printer
//        sp_tracker->SetPrinter( sp_printer);
//      }

      // create the temperature control
      util::ShPtr< mc::TemperatureInterface> sp_temperature
      (
        new mc::TemperatureAccepted
        (
          mc::TemperatureAccepted::GetParameterStartFraction()->GetNumericalValue< double>(),
          mc::TemperatureAccepted::GetParameterEndFraction()->GetNumericalValue< double>(),
          m_MaxIterations,
          mc::TemperatureAccepted::GetParameterStartTemperature()->GetNumericalValue< double>(),
          mc::TemperatureAccepted::GetParameterUpdateInterval()->GetNumericalValue< size_t>()
        )
      );

      // create the metropolis
      mc::Metropolis< double> metropolis( sp_temperature, true);

      // create the termination criterion
      opti::CriterionCombine< assemble::ProteinModel, double> termination_criterion;
      termination_criterion.InsertCriteria
      (
        opti::CriterionNumberIterations< assemble::ProteinModel, double>
        (
          m_MaxIterations
        )
      );
      termination_criterion.InsertCriteria
      (
        opti::CriterionUnimproved< assemble::ProteinModel, double>
        (
          m_MaxIterations / s_MaxIterationsStepsInARowFraction
        )
      );

      // create the approximator
      mc::Approximator< assemble::ProteinModel, double> approximator
      (
        *density_agreement,
        mutate,
        metropolis,
        termination_criterion,
        PROTEIN_MODEL
      );

      // approximate
      approximator.Approximate();

      return approximator.GetTracker().GetBest()->First();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FitProteinMinimizerMC::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FitProteinMinimizerMC::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace density
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
#include "density/bcl_density_fit_protein_minimizer_powell.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "density/bcl_density_map.h"
#include "opti/bcl_opti_approximator_powell.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_number_iterations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

  ///////////
  // types //
  ///////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FitProteinMinimizerPowell::PositionCorrelation::PositionCorrelation
    (
      const util::ShPtr< ProteinAgreementInterface> &AGREEMENT,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) :
      m_Agreement( AGREEMENT),
      m_Model( util::ToSiPtr( PROTEIN_MODEL))
    {
    }

    //! @brief Clone function
    //! @return pointer to new PositionCorrelation
    FitProteinMinimizerPowell::PositionCorrelation *FitProteinMinimizerPowell::PositionCorrelation::Clone() const
    {
      return new PositionCorrelation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FitProteinMinimizerPowell::PositionCorrelation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the correlation
    //! @param VECTOR containing 6 elements - three rotations, three translations
    double FitProteinMinimizerPowell::PositionCorrelation::operator()( const linal::Vector< double> &VECTOR) const
    {
      // return the agreement
      return m_Agreement->operator ()( *TransformedHardCopy( *m_Model, VECTOR));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FitProteinMinimizerPowell::PositionCorrelation::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FitProteinMinimizerPowell::PositionCorrelation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    util::ShPtr< assemble::ProteinModel> FitProteinMinimizerPowell::PositionCorrelation::TransformedHardCopy
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const linal::Vector< double> &VECTOR
    )
    {
      // copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.HardCopy());

      // create transformation
      const math::TransformationMatrix3D transformation( VECTOR);

      // apply transformation
      new_model->Transform( transformation);

      // end
      return new_model;
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FitProteinMinimizerPowell::s_Instance
    (
      GetObjectInstances().AddInstance( new FitProteinMinimizerPowell())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FitProteinMinimizerPowell::FitProteinMinimizerPowell()
    {
    }

    //! @brief Clone function
    //! @return pointer to new FitProteinMinimizerPowell
    FitProteinMinimizerPowell *FitProteinMinimizerPowell::Clone() const
    {
      return new FitProteinMinimizerPowell( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FitProteinMinimizerPowell::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the resolution of the density map
    //! @param RESOLUTION density map and simulation resolution
    void FitProteinMinimizerPowell::SetResolution( const double RESOLUTION)
    {
      m_Resolution = RESOLUTION;
    }

    //! @brief set max translation and rotation
    //! @param MAX_TRANSLATION max translation in any direction for a single iteration
    //! @param MAX_ROTATION max rotation in radians in any direction for a single iteration
    void FitProteinMinimizerPowell::SetMaxTranslationAndRotation( const double MAX_TRANSLATION, const double MAX_ROTATION)
    {
      m_MaxTranslation = MAX_TRANSLATION;
      m_MaxRotation    = MAX_ROTATION;
    }

    //! @brief set the max number of iterations for minimization
    //! @param MAX_NUMBER_ITERATIONS maximum number of iterations for minimization
    void FitProteinMinimizerPowell::SetMaxIterations( const size_t MAX_NUMBER_ITERATIONS)
    {
      m_MaxIterations = MAX_NUMBER_ITERATIONS;
    }

    //! @brief set protein agreement measure to be used
    //! @param AGREEMENT protein agreement enumerator
    void FitProteinMinimizerPowell::SetProteinAgreement( const ProteinAgreement &AGREEMENT)
    {
      m_ProteinAgreement = AGREEMENT;
    }

    //! @brief simulator to use
    //! @param DENSITY_SIMULATOR simulator enumerator
    void FitProteinMinimizerPowell::SetSimulator( const Simulator &DENSITY_SIMULATOR)
    {
      m_Simulator = DENSITY_SIMULATOR;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator minimizing the position of a protein model within a given density map
    //! @param PROTEIN_MODEL start position of given protein model
    //! @param DENSITY_MAP the density map to fit the PROTEIN_MODEL into
    //! @return the fitted protein model
    assemble::ProteinModel FitProteinMinimizerPowell::operator()( const assemble::ProteinModel &PROTEIN_MODEL, const Map &DENSITY_MAP) const
    {
      // construct a function that calculates the deviation between simulated density from atoms and the given density map
      util::ShPtr< ProteinAgreementInterface> density_agreement
      (
        GetProteinAgreements().CreateProteinAgreement
        (
          m_ProteinAgreement,
          m_Simulator,
          util::ToSiPtr( DENSITY_MAP),
          m_Resolution
        )
      );

      // create function opencl
      util::ShPtr< math::FunctionInterfaceSerializable< linal::Vector< double>, double> > sp_objective
      (
        new PositionCorrelation
        (
          density_agreement,
          PROTEIN_MODEL
        )
      );

      storage::Vector< linal::Vector< double> > search_directions;
      linal::Vector< double> direction( 6, 0.0);
      const linal::Vector< double> start( direction);
      direction( 0) = m_MaxRotation;
      search_directions.PushBack( direction);
      direction( 0) = 0.0;
      direction( 1) = m_MaxRotation;
      search_directions.PushBack( direction);
      direction( 1) = 0.0;
      direction( 2) = m_MaxRotation * 0.5;
      search_directions.PushBack( direction);
      direction( 2) = 0.0;
      direction( 3) = m_MaxTranslation;
      search_directions.PushBack( direction);
      direction( 3) = 0.0;
      direction( 4) = m_MaxTranslation;
      search_directions.PushBack( direction);
      direction( 4) = 0.0;
      direction( 5) = m_MaxTranslation;
      search_directions.PushBack( direction);

      // create termination criteria for the approximation
      opti::CriterionCombine< linal::Vector< double>, double> criterion_combine;
      criterion_combine.InsertCriteria
      (
        opti::CriterionNumberIterations< linal::Vector< double>, double>( m_MaxIterations / 10)
      );
      criterion_combine.InsertCriteria
      (
        opti::CriterionConvergenceResult< linal::Vector< double>, double>( 1, 0.001)
      );

      // create powell approximator from its members
      opti::ApproximatorPowell< linal::Vector< double>, double> approximator
      (
        *sp_objective, criterion_combine, search_directions, start
      );

      // do the actual approximation
      approximator.Approximate();

      // create the final model
      const util::ShPtr< assemble::ProteinModel> best_model
      (
        PositionCorrelation::TransformedHardCopy( PROTEIN_MODEL, approximator.GetTracker().GetBest()->First())
      );

      return *best_model;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FitProteinMinimizerPowell::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FitProteinMinimizerPowell::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace density
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
#include "density/bcl_density_fit_protein_minimizers.h"

// includes from bcl - sorted alphabetically
#include "density/bcl_density_fit_protein_minimizer_mc.h"
#include "density/bcl_density_fit_protein_minimizer_powell.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FitProteinMinimizers::FitProteinMinimizers() :
      e_MC(     AddEnum( "MC"    , util::ShPtr< FitProteinMinimizerInterface>( new FitProteinMinimizerMC()))),
      e_Powell( AddEnum( "Powell", util::ShPtr< FitProteinMinimizerInterface>( new FitProteinMinimizerPowell())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FitProteinMinimizers::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief construct on access function for all FitProteinMinimizers
    //! @return reference to only instances of FitProteinMinimizers
    FitProteinMinimizers &GetFitProteinMinimizers()
    {
      return FitProteinMinimizers::GetEnums();
    }

  } // namespace density

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< density::FitProteinMinimizerInterface>, density::FitProteinMinimizers>;

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
#include "density/bcl_density_map.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_point_cloud.h"
#include "math/bcl_math_angle.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Map::s_Instance
    (
      GetObjectInstances().AddInstance( new Map())
    );

    //! @brief default angle
    //! @return Vector3D with three default angles
    const linal::Vector3D &Map::GetDefaultAngle()
    {
      // default is 90 degrees
      static const linal::Vector3D s_angle( 90.0, 90.0, 90.0);
      return s_angle;
    }

    //! @brief default axis
    //! @brief storage vector ND indicating indes for slow, middle and fast changing axis
    const linal::VectorND< int, 3> &Map::GetDefaultAxis()
    {
      // default is first index for slow changing index, 3 for fast chaning index
      static const linal::VectorND< int, 3> s_axis( 1, 2, 3);
      return s_axis;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Map::Map() :
      m_Index(),
      m_Intervals(),
      m_Length(),
      m_CellWidth(),
      m_Angle( GetDefaultAngle()),
      m_Axis( GetDefaultAxis()),
      m_Minimum(),
      m_Maximum(),
      m_Mean(),
      m_SpaceGroup( 0),
      m_NumberBytesSymmetryData( 0),
      m_Rmsd(),
      m_Origin(),
      m_Data(),
      m_MachineStamp( 0)
    {
      std::fill_n( m_Extra, s_ExtraSize, '0');
    }

    //! @brief construct from all map parameters
    //! @brief DATA the data as tensor
    Map::Map
    (
      const math::Tensor< double> &DATA,
      const linal::VectorND< int, 3> &INDEX,
      const linal::VectorND< int, 3> &INTERVALS,
      const linal::Vector3D &LENGTH,
      const linal::Vector3D &CELLWIDTH,
      const linal::Vector3D &ANGLE,
      const linal::VectorND< int, 3> &AXIS,
      const linal::Vector3D &ORIGIN
    ) :
      m_Index( INDEX),
      m_Intervals( INTERVALS),
      m_Length( LENGTH),
      m_CellWidth( CELLWIDTH),
      m_Angle( ANGLE),
      m_Axis( AXIS),
      m_Minimum(),
      m_Maximum(),
      m_Mean(),
      m_SpaceGroup( 0),
      m_NumberBytesSymmetryData( 0),
      m_Rmsd(),
      m_Origin( ORIGIN),
      m_Data( DATA),
      m_MachineStamp( 0),
      m_Labels()
    {
      std::fill_n( m_Extra, s_ExtraSize, '0');
      CalculateMinMaxMeanRmsd();
    }

    //! @brief copy constructor
    //! @param DENSITY_MAP map to copy from
    Map::Map( const Map &DENSITY_MAP) :
      m_Index( DENSITY_MAP.m_Index),
      m_Intervals( DENSITY_MAP.m_Intervals),
      m_Length( DENSITY_MAP.m_Length),
      m_CellWidth( DENSITY_MAP.m_CellWidth),
      m_Angle( DENSITY_MAP.m_Angle),
      m_Axis( DENSITY_MAP.m_Axis),
      m_Minimum( DENSITY_MAP.m_Minimum),
      m_Maximum( DENSITY_MAP.m_Maximum),
      m_Mean( DENSITY_MAP.m_Mean),
      m_SpaceGroup( DENSITY_MAP.m_SpaceGroup),
      m_NumberBytesSymmetryData( DENSITY_MAP.m_NumberBytesSymmetryData),
      m_Rmsd( DENSITY_MAP.m_Rmsd),
      m_Origin( DENSITY_MAP.m_Origin),
      m_Data( DENSITY_MAP.m_Data),
      m_MachineStamp( DENSITY_MAP.m_MachineStamp),
      m_Labels( DENSITY_MAP.m_Labels)
    {
      std::copy( DENSITY_MAP.m_Extra, DENSITY_MAP.m_Extra + s_ExtraSize, m_Extra);
    }

  ////////////////
  // operations //
  ////////////////

    //! cut out a part of the density map
    //! @param POSCOL   column or x position
    //! @param POSROW   row    pr y position
    //! @param POSLAYER layer  or z position
    //! @param NCOL     number of columns
    //! @param NROW     number of rows
    //! @param NLAYER   number of layers
    //! @return sub density map from given position of given size
    Map Map::SubMap
    (
      const size_t POSCOL, const size_t POSROW, const size_t POSLAYER,
      const size_t NCOL, const size_t NROW, const size_t NLAYER
    ) const
    {
      //make sure that coordinates does not extent the Denistymap
      BCL_Assert
      (
           POSCOL   < m_Data.GetNumberCols()   && POSCOL   + NCOL   <= m_Data.GetNumberCols()
        && POSROW   < m_Data.GetNumberRows()   && POSROW   + NROW   <= m_Data.GetNumberRows()
        && POSLAYER < m_Data.NumberLayers() && POSLAYER + NLAYER <= m_Data.NumberLayers(),
        "Sub coordinates extent Densitymap"
      );

      // extension output map
      const storage::VectorND< 3, size_t> ext_output( NCOL, NROW, NLAYER);

      // index changes
      linal::VectorND< int, 3> new_index( 0, 0, 0);
      new_index( 0) = m_Index( 0) + POSCOL;
      new_index( 1) = m_Index( 1) + POSROW;
      new_index( 2) = m_Index( 2) + POSLAYER;

      // intervals
      linal::VectorND< int, 3> intervals( 0, 0, 0);
      for( size_t i( 0); i < 3; ++i)
      {
        if( new_index( i) <= -int( ext_output( i)))
        {
          intervals( i) = math::Absolute( new_index( i));
        }
        else if( new_index( i) >= 0)
        {
          intervals( i) = new_index( i) + ext_output( i) - 1;
        }
        else
        {
          intervals( i) = ext_output( i) - 1;
        }
      }

      // length changes
      linal::Vector3D new_length( 0.0);
      new_length( 0) = intervals( 0) * m_CellWidth( 0);
      new_length( 1) = intervals( 1) * m_CellWidth( 1);
      new_length( 2) = intervals( 2) * m_CellWidth( 2);

      // create submap and return
      return Map
      (
        m_Data.SubTensor( NLAYER, NROW, NCOL, POSLAYER, POSROW, POSCOL),
        new_index,
        intervals,
        new_length,
        m_CellWidth,
        m_Angle,
        m_Axis,
        m_Origin
      );
    }

    //! get the common sub tensor from this and argument density map
    //! @param DENSITY_MAP the second density map to find common subtensor
    //! @return two tensors, first tensor of this map, second tensot of argument
    storage::VectorND< 2, math::Tensor< double> > Map::CommonSubTensor( const Map &DENSITY_MAP) const
    {
      // make sure that the two density maps have the same CellWidth
      BCL_Assert
      (
        math::EqualWithinTolerance( m_CellWidth, DENSITY_MAP.m_CellWidth),
        "Calculating common sub tensor requires density maps with equal CellWidth\nthis: " +
        util::Format()( m_CellWidth) + "\nargument: " + util::Format()( DENSITY_MAP.m_CellWidth)
      );

      // make sure that the origins are the same
      BCL_Assert
      (
        math::EqualWithinTolerance( m_Origin, DENSITY_MAP.m_Origin),
        "Calculating common sub tensor requires density maps with equal Origin\nthis: " +
        util::Format()( m_Origin) + "\nargument: " + util::Format()( DENSITY_MAP.m_Origin)
      );

      // relative index of argument density relative to this
      linal::Vector< int> this_start( 3, int( 0));
      linal::Vector< int> arg_start(  3, int( 0));
      linal::Vector< int> this_end(   3, int( 0));
      linal::Vector< int> arg_end(    3, int( 0));
      // find common start and end
      linal::Vector< int> common_start( 3, int( 0));
      linal::Vector< int> common_end( 3, int( 0));
      linal::Vector< int> extent( 3, int( 0));

      // iterate over dimensions
      for( size_t i( 0); i < 3; ++i)
      {
        this_start( i) =             m_Index( i);
        arg_start(  i) = DENSITY_MAP.m_Index( i);
        const int this_end      = this_start( i) +             GetDimensions()( i);
        const int arg_end       = arg_start(  i) + DENSITY_MAP.GetDimensions()( i);
        const int common_start  = std::max( this_start( i), arg_start( i));
        const int common_end    = std::min( this_end      , arg_end);
        extent( i) = common_end - common_start - 1;

        // no common sub density
        if( extent( i) <= 0)
        {
          return storage::VectorND< 2, math::Tensor< double> >();
        }

        this_start( i) = common_start -             m_Index( i);
        arg_start(  i) = common_start - DENSITY_MAP.m_Index( i);
      }

      // create sub tensors and return
      return storage::VectorND< 2, math::Tensor< double> >
      (
        m_Data.SubTensor( extent( 2), extent( 1), extent( 0), this_start( 2), this_start( 1), this_start( 0)),
        DENSITY_MAP.m_Data.SubTensor( extent( 2), extent( 1), extent( 0), arg_start( 2), arg_start( 1), arg_start( 0))
      );
    }

    //! @brief calculate the min max and mean from the intensities stored in the data
    void Map::CalculateMinMaxMeanRmsd()
    {
      m_Minimum = math::Statistics::MinimumValue( m_Data.Begin(), m_Data.End());
      m_Maximum = math::Statistics::MaximumValue( m_Data.Begin(), m_Data.End());
      m_Mean    = math::Statistics::Mean( m_Data.Begin(), m_Data.End());
      m_Rmsd    = math::Statistics::StandardDeviation( m_Data.Begin(), m_Data.End());
    }

    //! compute histogram over all densities in BIN bins
    math::Histogram Map::Histogram( const size_t NR_BINS) const
    {
      math::Histogram histogram( m_Minimum, ( m_Maximum - m_Minimum) / NR_BINS, NR_BINS);
      for( const double *ptr( m_Data.Begin()), *ptr_end( m_Data.End()); ptr != ptr_end; ++ptr)
      {
        histogram.PushBack( *ptr);
      }

      return histogram;
    }

    //! compute a point cloud of NUMBER_OF_POINTS with an at least distance of FEATURE_DISTANCE in Angstrom using
    //! the highest derivatives between voxels in the map and the own intensity as the user can weight
    coord::PointCloud
    Map::CalculatePointCloud
    (
      const size_t NUMBER_OF_POINTS,
      const double FEATURE_DISTANCE,
      const double RATIO_INTENSITY_GRADIENT
    ) const
    {
      // local variables
      // function, that evaluates if a voxel overlaps by the given radius with another voxel
      const VoxelOverlap does_voxel_overlap_within_radius
      (
        FEATURE_DISTANCE,
        m_CellWidth
      );

      // set that keeps track of voxels with highest intensities and gradients
      std::multiset< Voxel, VoxelDensityGreater> considered_voxels;

      // set that stores the voxels that do overlap, but might be not considered too early
      std::multiset< Voxel, VoxelDensityGreater> overlapping_voxels;

      // store min intensity
      const double min_intensity( GetMean());

      // min intensity of a voxels has to be at least 50% of the highest intensity
      BCL_MessageVrb( "PointCloud( NUMBER_OF_POINTS = " + util::Format()( NUMBER_OF_POINTS) +
        ", FEATURE_DISTANCE = " + util::Format()( FEATURE_DISTANCE) + ")\nstart of building Pointcloud");

      // loop over all voxels and memorize best NUMBER_OF_POINTS
      for( size_t i( 1); i < m_Data.NumberLayers() - 1; ++i)
      {
        for( size_t j( 1); j < m_Data.GetNumberRows() - 1; ++j)
        {
          for( size_t k( 1); k < m_Data.GetNumberCols() - 1; ++k)
          {
            const double current_intensity( m_Data( i, j, k));

            if( current_intensity < min_intensity)
            {
              continue;
            }

            // make new voxel
            const double derivative_intensity
            (
              RATIO_INTENSITY_GRADIENT * current_intensity +
              double( 1) / double( 24) *
              (
                // direct nighbors
                math::Absolute( current_intensity - m_Data( i - 1, j    , k    )) +
                math::Absolute( current_intensity - m_Data( i + 1, j    , k    )) +
                math::Absolute( current_intensity - m_Data( i    , j - 1, k    )) +
                math::Absolute( current_intensity - m_Data( i    , j + 1, k    )) +
                math::Absolute( current_intensity - m_Data( i    , j    , k - 1)) +
                math::Absolute( current_intensity - m_Data( i    , j    , k + 1))
              ) +
              double( 1) / ( double( 24) * math::Sqrt( double( 2))) *
              (
                // neighbors on edges
                math::Absolute( current_intensity - m_Data( i - 1, j - 1, k    )) +
                math::Absolute( current_intensity - m_Data( i - 1, j + 1, k    )) +
                math::Absolute( current_intensity - m_Data( i + 1, j - 1, k    )) +
                math::Absolute( current_intensity - m_Data( i + 1, j + 1, k    )) +
                math::Absolute( current_intensity - m_Data( i    , j - 1, k - 1)) +
                math::Absolute( current_intensity - m_Data( i    , j - 1, k + 1)) +
                math::Absolute( current_intensity - m_Data( i    , j + 1, k - 1)) +
                math::Absolute( current_intensity - m_Data( i    , j + 1, k + 1)) +
                math::Absolute( current_intensity - m_Data( i - 1, j    , k - 1)) +
                math::Absolute( current_intensity - m_Data( i - 1, j    , k + 1)) +
                math::Absolute( current_intensity - m_Data( i + 1, j    , k - 1)) +
                math::Absolute( current_intensity - m_Data( i + 1, j    , k + 1))
              ) +
              double( 1) / ( double( 24) * math::Sqrt( double( 3))) *
              (
                // neighbor on corners
                math::Absolute( current_intensity - m_Data( i - 1, j - 1, k - 1)) +
                math::Absolute( current_intensity - m_Data( i - 1, j - 1, k + 1)) +
                math::Absolute( current_intensity - m_Data( i - 1, j + 1, k - 1)) +
                math::Absolute( current_intensity - m_Data( i - 1, j + 1, k + 1)) +
                math::Absolute( current_intensity - m_Data( i + 1, j - 1, k - 1)) +
                math::Absolute( current_intensity - m_Data( i + 1, j - 1, k + 1)) +
                math::Absolute( current_intensity - m_Data( i + 1, j + 1, k - 1)) +
                math::Absolute( current_intensity - m_Data( i + 1, j + 1, k + 1))
              )
            );
            Voxel current_voxel = { i, j, k, derivative_intensity};

            std::vector< std::multiset< Voxel, VoxelDensityGreater>::iterator> current_overlaps;
            // find all overlapping voxels in the set of considered voxels
            std::multiset< Voxel, VoxelDensityGreater>::iterator
              voxel_itr_overlap( considered_voxels.begin()), voxel_itr_end( considered_voxels.end());
            bool all_smaller_intensity( true);
            while( true)
            {
              voxel_itr_overlap =
                std::find_if
                (
                  voxel_itr_overlap,
                  voxel_itr_end,
                  std::bind2nd( does_voxel_overlap_within_radius, current_voxel)
                );
              if( voxel_itr_overlap == voxel_itr_end)
              {
                break;
              }
              else
              {
                current_overlaps.push_back( voxel_itr_overlap);
                all_smaller_intensity &= voxel_itr_overlap->m_Intensity < current_voxel.m_Intensity;
              }
              ++voxel_itr_overlap;
            }

            // insert if no overlap was found
            if( current_overlaps.empty())
            {
              considered_voxels.insert( current_voxel);
            }
            else
            {
              // replace with the higher intensity
              if( all_smaller_intensity)
              {
                considered_voxels.insert( current_voxel);
                for
                (
                  std::vector< std::multiset< Voxel, VoxelDensityGreater>::iterator>::const_iterator
                    cur_overlap_itr( current_overlaps.begin()), cur_overlap_itr_end( current_overlaps.end());
                  cur_overlap_itr != cur_overlap_itr_end;
                  ++cur_overlap_itr
                )
                {
                  overlapping_voxels.insert( **cur_overlap_itr);
                  considered_voxels.erase( *cur_overlap_itr);
                }
              }
              // insert it in the overlapping voxels
              else
              {
                overlapping_voxels.insert( current_voxel);
              }
            }
            // check if considered voxels exceeds the considered size
            if( considered_voxels.size() > NUMBER_OF_POINTS)
            {
              std::multiset< Voxel, VoxelDensityGreater>::iterator itr( considered_voxels.end());
              --itr;
              overlapping_voxels.insert( *itr);
              considered_voxels.erase( itr);
            }
          }
          // remove overlapping voxels that will not be used any more
          if( !overlapping_voxels.empty() && !considered_voxels.empty())
          {
            // find first voxel in overlapping voxels that has higher intensity than the worst in the considered set
            std::multiset< Voxel, VoxelDensityGreater>::iterator voxel_itr_end( overlapping_voxels.end()), voxel_itr_overlap
            (
              std::find_if( overlapping_voxels.begin(), voxel_itr_end, std::bind1st( VoxelDensityGreater(), *( --considered_voxels.end())))
            );
            // erase the range with intensity that is too low
            overlapping_voxels.erase( voxel_itr_overlap, voxel_itr_end);
          }
        }
      }

      // if an overlapping Voxel has a higher intensity than one of the considered Voxels and non of the Voxels that
      // have a higher intensity than this, it has to be inserted and all overlapping Voxels with lower intensity have
      // to be erased

      // iterator on the overlapping voxels
      std::multiset< Voxel, VoxelDensityGreater>::iterator
        over_voxel_itr( overlapping_voxels.begin()), over_voxel_itr_end( overlapping_voxels.end());

      // iterator on the considered
      std::multiset< Voxel, VoxelDensityGreater>::iterator
        cons_voxel_itr_low( considered_voxels.begin());

      // iterate as long as there are overlapping voxels
      while( over_voxel_itr != over_voxel_itr_end)
      {
        std::multiset< Voxel, VoxelDensityGreater>::iterator cons_voxel_itr_end( considered_voxels.end());
        // find the first voxel with lower intensity in the considered voxels
        std::multiset< Voxel, VoxelDensityGreater>::iterator cons_voxel_itr
        (
          std::find_if( cons_voxel_itr_low, cons_voxel_itr_end, std::bind1st( VoxelDensityGreater(), *over_voxel_itr))
        );

        // no voxel with lower intensity
        if( cons_voxel_itr == cons_voxel_itr_end)
        {
          break;
        }

        // search in the range that has a higher intensity for an overlapping voxel
        std::multiset< Voxel, VoxelDensityGreater>::iterator cons_higher_overlapping_voxel_itr
        (
          std::find_if( considered_voxels.begin(), cons_voxel_itr, std::bind2nd( does_voxel_overlap_within_radius, *over_voxel_itr))
        );

        // if there is no overlapping iterator in that range
        if( cons_higher_overlapping_voxel_itr == cons_voxel_itr)
        {
          // insert that voxel
          considered_voxels.insert( *over_voxel_itr);
        }
        // proceed to next
        ++over_voxel_itr;
      }

      // delete all overlapping voxels with lower intensities
      std::multiset< Voxel, VoxelDensityGreater>::iterator voxel_itr( considered_voxels.begin()), voxel_itr_end( considered_voxels.end());
      while( voxel_itr != voxel_itr_end)
      {
        std::multiset< Voxel, VoxelDensityGreater>::iterator voxel_itr_next( voxel_itr);
        ++voxel_itr_next;
        // at end
        if( voxel_itr_next == voxel_itr_end)
        {
          break;
        }

        // delete all overlapping points
        while( true)
        {
          std::multiset< Voxel, VoxelDensityGreater>::iterator voxel_itr_overlap
          (
            std::find_if( voxel_itr_next, voxel_itr_end, std::bind2nd( does_voxel_overlap_within_radius, *voxel_itr))
          );

          // remove overlapping
          if( voxel_itr_overlap != voxel_itr_end)
          {
            // progress after the overlapping iterator
            voxel_itr_next = voxel_itr_overlap;
            ++voxel_itr_next;
            considered_voxels.erase( voxel_itr_overlap);
          }
          // no remaining overlapping voxel
          else
          {
            break;
          }
        }
        ++voxel_itr;
      }

      BCL_MessageVrb
      (
        "NUMBER_OF_POINTS = " + util::Format()( NUMBER_OF_POINTS) + ", FEATURE_DISTANCE = " + util::Format()( FEATURE_DISTANCE) +
        "\nPoints found at all: " + util::Format()( considered_voxels.size()) + " will be reduced to target number of points"
      );

      // build PointCloud
      coord::PointCloud pointcloud;
      pointcloud.GetData().AllocateMemory( NUMBER_OF_POINTS);

      size_t count( 0);
      for
      (
        std::multiset< Voxel, VoxelDensityGreater>::const_iterator
          itr( considered_voxels.begin()), itr_end( considered_voxels.end());
        itr != itr_end && count < NUMBER_OF_POINTS;
        ++itr, ++count
      )
      {
        pointcloud.PushBack
        (
          linal::Vector3D
          (
            double( ( int( itr->m_Col)   + m_Index( 0)) - 0.5) * m_CellWidth.X() + m_Origin.X(),
            double( ( int( itr->m_Row)   + m_Index( 1)) - 0.5) * m_CellWidth.Y() + m_Origin.Y(),
            double( ( int( itr->m_Layer) + m_Index( 2)) - 0.5) * m_CellWidth.Z() + m_Origin.Z()
          )
        );
      }

      return pointcloud;
    }

    //! @brief add noise
    //! @param RNG an random number Generator object
    //! @param MEAN the mean of the gaussian distribution
    //! @param STANDARD_DEVIATION standard deviation of the gaussian distribution
    //! @return the cross correlation coefficient map (current map as simulated) vs new map with noise
    double Map::AddNoise( const random::DistributionInterface &RNG, const double MEAN, const double STANDARD_DEVIATION)
    {
      // copy of this map for cross correlation calculation after noise addition
      Map original_map( *this);

      // iterate over map and add random noise
      for( double *ptr( m_Data.Begin()), *ptr_end( m_Data.End()); ptr != ptr_end; ++ptr)
      {
        *ptr += RNG.RandomGaussian( MEAN, STANDARD_DEVIATION);
      }

      // recalculate min, max, mean and rmsd
      CalculateMinMaxMeanRmsd();

      // end - calculate the correlation to the original map
      return CrossCorrelationCoefficient( original_map, original_map.m_Minimum);
    }

    //! @brief normalize
    //! transforms intensities, so that rmsd is 1
    void Map::Normalize()
    {
      // subtract mean from all intensities and devide by rmsd
      m_Data /= m_Rmsd;

      // recalculate min, max, mean and rmsd
      CalculateMinMaxMeanRmsd();
    }

    //! calculate the standard deviation to a smaller density map
    //! http://en.wikipedia.org/wiki/Cross-correlation#Normalized_cross-correlation
    double Map::Correlation( const Map &DENSITY_MAP) const
    {
      return CrossCorrelationCoefficient( DENSITY_MAP, double( 0));
    }

    //! calculate the cross correlation factor to a simulated density map, only using voxels that are above the given
    //! CONTOUR_LEVEL in the simulated density map
    //! @param SIMULATED_DENSITY_MAP
    //! @param CONTOUR_LEVEL
    //! @return cross correlation coefficient
    double Map::CrossCorrelationCoefficient
    (
      const Map &SIMULATED_DENSITY_MAP,
      const double CONTOUR_LEVEL
    ) const
    {
      // make sure that the two density maps have the same CellWidth
      BCL_Assert
      (
        math::EqualWithinTolerance( m_CellWidth, SIMULATED_DENSITY_MAP.m_CellWidth),
        "Calculating cross correlation factor requires density maps with equal CellWidth\nthis: " +
        util::Format()( m_CellWidth) + "\nargument: " + util::Format()( SIMULATED_DENSITY_MAP.m_CellWidth)
      );

      //store the index for the lower right and upper left corner of the argument density map
      linal::Vector< int> lower_right_corner_argument( 3, int( 0)), upper_left_corner_argument( 3, int( 0));

      for( size_t i( 0); i < 3; ++i)
      {
        lower_right_corner_argument( i) = int( ( SIMULATED_DENSITY_MAP.m_Origin( i) - m_Origin( i)) / m_CellWidth( i) + double( SIMULATED_DENSITY_MAP.m_Index( i) - m_Index( i)));
        upper_left_corner_argument( i) = lower_right_corner_argument( i) + SIMULATED_DENSITY_MAP.GetDimensions()( i) - 1;
      }

      // make sure that at all indices of any corner of the argument density map is within the Map
      if(
          !(
                upper_left_corner_argument(  0) >= 0 && lower_right_corner_argument( 0) <= int( GetDimensions()( 0))
             && upper_left_corner_argument(  1) >= 0 && lower_right_corner_argument( 1) <= int( GetDimensions()( 1))
             && upper_left_corner_argument(  2) >= 0 && lower_right_corner_argument( 2) <= int( GetDimensions()( 2))
            )
        )
      {
        return util::GetUndefined< double>();
      }

      // number of voxels that are above the CONTOUR_LEVEL in the experimental (this) density map
      size_t count_voxel( 0);

      // mean and standard deviation
      math::RunningAverageSD< double> mean_sd_this;
      math::RunningAverageSD< double> mean_sd_sim;

      // define the indices to iterate over
      const size_t i1_min( size_t( std::max( int( 0),  lower_right_corner_argument( 0))));
      const size_t i2_min( size_t( std::max( int( 0), -lower_right_corner_argument( 0))));
      const size_t i1_max( GetDimensions()( 0));
      const size_t i2_max( SIMULATED_DENSITY_MAP.GetDimensions()( 0));

      const size_t j1_min( size_t( std::max( int( 0),  lower_right_corner_argument( 1))));
      const size_t j2_min( size_t( std::max( int( 0), -lower_right_corner_argument( 1))));
      const size_t j1_max( GetDimensions()( 1));
      const size_t j2_max( SIMULATED_DENSITY_MAP.GetDimensions()( 1));

      const size_t k1_min( size_t( std::max( int( 0),  lower_right_corner_argument( 2))));
      const size_t k2_min( size_t( std::max( int( 0), -lower_right_corner_argument( 2))));
      const size_t k1_max( GetDimensions()( 2));
      const size_t k2_max( SIMULATED_DENSITY_MAP.GetDimensions()( 2));

      // calculate correlation by iterating over all density values
      for( size_t i1( i1_min), i2( i2_min); i1 < i1_max && i2 < i2_max; ++i1, ++i2)
      {
        for( size_t j1( j1_min), j2( j2_min); j1 < j1_max && j2 < j2_max; ++j1, ++j2)
        {
          for( size_t k1( k1_min), k2( k2_min); k1 < k1_max && k2 < k2_max; ++k1, ++k2)
          {
            const double this_int( operator()( i1, j1, k1));

            // skip exp intensities below contour level
            if( this_int < CONTOUR_LEVEL)
            {
              continue;
            }

            const double sim_int( SIMULATED_DENSITY_MAP( i2, j2, k2));
            ++count_voxel;
            mean_sd_this += this_int;
            mean_sd_sim += sim_int;
          }
        }
      }

      const double mean_this( mean_sd_this.GetAverage());
      const double mean_sim( mean_sd_sim.GetAverage());

      double correlation( 0);

      // calculate correlation by iterating over all density values
      for( size_t i1( i1_min), i2( i2_min); i1 < i1_max && i2 < i2_max; ++i1, ++i2)
      {
        for( size_t j1( j1_min), j2( j2_min); j1 < j1_max && j2 < j2_max; ++j1, ++j2)
        {
          for( size_t k1( k1_min), k2( k2_min); k1 < k1_max && k2 < k2_max; ++k1, ++k2)
          {
            const double this_int( operator()( i1, j1, k1));

            // skip exp intensities below contour level
            if( this_int < CONTOUR_LEVEL)
            {
              continue;
            }

            const double sim_int( SIMULATED_DENSITY_MAP( i2, j2, k2));
            correlation += ( this_int - mean_this) * ( sim_int - mean_sim);
          }
        }
      }

      // normalize by number of voxel
      correlation = double( 1) / ( double( count_voxel * mean_sd_this.GetStandardDeviation() * mean_sd_sim.GetStandardDeviation())) * correlation;

      return correlation;
    }

    //! calculate the standard deviation to a smaller density map
    double Map::StandardDeviation( const Map &DENSITY_MAP) const
    {
      // make sure that the two density maps have the same CellWidth
      BCL_Assert
      (
        m_CellWidth == DENSITY_MAP.m_CellWidth,
        "Calculating correlation factor requires density maps with equal CellWidth"
      );

      //store the index for the lower right and upper left corner of the argument density map
      linal::Vector< int> lower_right_corner_argument( 3, int( 0)), upper_left_corner_argument( 3, int( 0));
      for( size_t i( 0); i < 3; ++i)
      {
        lower_right_corner_argument( i) = int( ( DENSITY_MAP.m_Origin( i) - m_Origin( i)) / m_CellWidth( i) + double( int( DENSITY_MAP.m_Index( i)) - int( m_Index( i))) + 0.5);
        upper_left_corner_argument( i) = lower_right_corner_argument( i) + DENSITY_MAP.GetDimensions()( i) - 1;
      }

      // make sure that at least one index of any corner of the argument density map is within the Map
      if(
          !(
                upper_left_corner_argument(  0) >= 0 && lower_right_corner_argument( 0) <= int( GetDimensions()( 0))
             && upper_left_corner_argument(  1) >= 0 && lower_right_corner_argument( 1) <= int( GetDimensions()( 1))
             && upper_left_corner_argument(  2) >= 0 && lower_right_corner_argument( 2) <= int( GetDimensions()( 2))
            )
        )
      {
        return util::GetUndefined< double>();
      }

      std::pair< double, double> sim_min_max_threshold_intensity( DENSITY_MAP.GetMinimum(), DENSITY_MAP.GetMaximum());

      // variables to store average value and standard deviation of experimanetal and model map
      // but only for occupied voxels in simulatet map
      std::pair< double, double> exp_avg_sd( 0, 0), sim_avg_sd( 0, 0);
      size_t count_occupied_voxels( 0), count_voxels( 0);

      //calculate average density over all occupied voxels in simulated map
      for( size_t i1( size_t( std::max( int( 0), lower_right_corner_argument( 2)))),
                  i2( size_t( std::max( int( 0), -lower_right_corner_argument( 2))));
                  i1 < GetDimensions()( 2) && i2 < DENSITY_MAP.GetDimensions()( 2); ++i1, ++i2)
        for( size_t j1( size_t( std::max( int( 0), lower_right_corner_argument( 1)))),
                    j2( size_t( std::max( int( 0), -lower_right_corner_argument( 1))));
                    j1 < GetDimensions()( 1) && j2 < DENSITY_MAP.GetDimensions()( 1); ++j1, ++j2)
          for( size_t k1( size_t( std::max( int( 0), lower_right_corner_argument( 0)))),
                      k2( size_t( std::max( int( 0), -lower_right_corner_argument( 0))));
                      k1 < GetDimensions()( 0) && k2 < DENSITY_MAP.GetDimensions()( 0); ++k1, ++k2)
            {
              exp_avg_sd.first += operator()( i1, j1, k1);
              sim_avg_sd.first += DENSITY_MAP( i2, j2, k2);
              count_voxels++;
            }
      exp_avg_sd.first /= count_voxels;
      sim_avg_sd.first /= count_voxels;

      BCL_MessageVrb( "calculating standard deviation for " + util::Format()( count_voxels) + " overlapping voxels");

      //calculate standard deviation in occupied voxels
      for( size_t i1( size_t( std::max( int( 0), lower_right_corner_argument( 2)))),
                  i2( size_t( std::max( int( 0), -lower_right_corner_argument( 2))));
                  i1 < GetDimensions()( 2) && i2 < DENSITY_MAP.GetDimensions()( 2); ++i1, ++i2)
        for( size_t j1( size_t( std::max( int( 0), lower_right_corner_argument( 1)))),
                    j2( size_t( std::max( int( 0), -lower_right_corner_argument( 1))));
                    j1 < GetDimensions()( 1) && j2 < DENSITY_MAP.GetDimensions()( 1); ++j1, ++j2)
          for( size_t k1( size_t( std::max( int( 0), lower_right_corner_argument( 0)))),
                      k2( size_t( std::max( int( 0), -lower_right_corner_argument( 0))));
                      k1 < GetDimensions()( 0) && k2 < DENSITY_MAP.GetDimensions()( 0); ++k1, ++k2)
            {
              exp_avg_sd.second += math::Sqr( exp_avg_sd.first - operator()( i1, j1, k1));
              sim_avg_sd.second += math::Sqr( sim_avg_sd.first - DENSITY_MAP( i2, j2, k2));
            }
      //calculate average density in Map and Argument
      exp_avg_sd.second /= ( count_voxels - 1);
      sim_avg_sd.second /= ( count_voxels - 1);

      //calculate standard deviation in Map and Argument
      exp_avg_sd.first = math::Sqrt( exp_avg_sd.first);
      sim_avg_sd.first = math::Sqrt( sim_avg_sd.first);
      sim_min_max_threshold_intensity.first  = sim_avg_sd.first - 3 * sim_avg_sd.second;
      sim_min_max_threshold_intensity.second = sim_avg_sd.first + 3 * sim_avg_sd.second;

      double standard_deviation( 0);

      // calculate standard deviation
      for( size_t i1( size_t( std::max( int( 0), lower_right_corner_argument( 2)))),
                  i2( size_t( std::max( int( 0), -lower_right_corner_argument( 2))));
                  i1 < GetDimensions()( 2) && i2 < DENSITY_MAP.GetDimensions()( 2); ++i1, ++i2)
      {
        for( size_t j1( size_t( std::max( int( 0), lower_right_corner_argument( 1)))),
                    j2( size_t( std::max( int( 0), -lower_right_corner_argument( 1))));
                    j1 < GetDimensions()( 1) && j2 < DENSITY_MAP.GetDimensions()( 1); ++j1, ++j2)
        {
          for( size_t k1( size_t( std::max( int( 0), lower_right_corner_argument( 0)))),
                      k2( size_t( std::max( int( 0), -lower_right_corner_argument( 0))));
                      k1 < GetDimensions()( 0) && k2 < DENSITY_MAP.GetDimensions()( 0); ++k1, ++k2)
          {
            // shift gaussian distributions to same avarage value and scale to same standard deviation and substract values in voxels but only if target density map voxel is occupied
//            if( DENSITY_MAP( i2, j2, k2) > sim_min_max_threshold_intensity.first && DENSITY_MAP( i2, j2, k2) < sim_min_max_threshold_intensity.second)
            {
              standard_deviation += math::Sqr( ( ( DENSITY_MAP( i2, j2, k2) - sim_avg_sd.first) * exp_avg_sd.second) / sim_avg_sd.second + exp_avg_sd.first - operator()( i1, j1, k1));
              count_occupied_voxels++;
            }
          }
        }
      }

      standard_deviation = math::Sqrt( standard_deviation / count_occupied_voxels);

      return standard_deviation;
    }

    //! @brief determine edges using a Sobel operator
    //! @sa http://en.wikipedia.org/wiki/Edge_detection
    Map &Map::EdgeDetectionSobel( const double WEIGTH_FACE, const double WEIGHT_EDGE, const double WEIGHT_VERTEX)
    {
      math::Tensor< double> edges( m_Data);
      edges.SetZero();

      const double s_mask_x[ 27] =
      {
         WEIGHT_VERTEX,  WEIGHT_EDGE,  WEIGHT_VERTEX,
         WEIGHT_EDGE,    WEIGTH_FACE,  WEIGHT_EDGE,
         WEIGHT_VERTEX,  WEIGHT_EDGE,  WEIGHT_VERTEX,

         0,  0,  0,
         0,  0,  0,
         0,  0,  0,

        -WEIGHT_VERTEX, -WEIGHT_EDGE, -WEIGHT_VERTEX,
        -WEIGHT_EDGE,   -WEIGTH_FACE, -WEIGHT_EDGE,
        -WEIGHT_VERTEX, -WEIGHT_EDGE, -WEIGHT_VERTEX
      };

      const double s_mask_y[ 27] =
      {
         WEIGHT_VERTEX,  0, -WEIGHT_VERTEX,
         WEIGHT_EDGE,    0, -WEIGHT_EDGE,
         WEIGHT_VERTEX,  0, -WEIGHT_VERTEX,

         WEIGHT_EDGE,  0, -WEIGHT_EDGE,
         WEIGTH_FACE,  0, -WEIGTH_FACE,
         WEIGHT_EDGE,  0, -WEIGHT_EDGE,

         WEIGHT_VERTEX,  0, -WEIGHT_VERTEX,
         WEIGHT_EDGE,    0, -WEIGHT_EDGE,
         WEIGHT_VERTEX,  0, -WEIGHT_VERTEX
      };

      const double s_mask_z[ 27] =
      {
         WEIGHT_VERTEX,  WEIGHT_EDGE,  WEIGHT_VERTEX,
         0, 0, 0,
        -WEIGHT_VERTEX, -WEIGHT_EDGE, -WEIGHT_VERTEX,

         WEIGHT_EDGE,  WEIGTH_FACE,  WEIGHT_EDGE,
         0, 0, 0,
        -WEIGHT_EDGE, -WEIGTH_FACE, -WEIGHT_EDGE,

         WEIGHT_VERTEX,  WEIGHT_EDGE,  WEIGHT_VERTEX,
         0, 0, 0,
        -WEIGHT_VERTEX, -WEIGHT_EDGE, -WEIGHT_VERTEX
      };

      const math::Tensor< double> mask_x( 3, 3, 3, s_mask_x);
      const math::Tensor< double> mask_y( 3, 3, 3, s_mask_y);
      const math::Tensor< double> mask_z( 3, 3, 3, s_mask_z);

      // iterate over all voxel
      for( size_t x = 1; x < m_Data.NumberLayers() - 1; ++x)
      {
        for( size_t y = 1; y < m_Data.GetNumberRows() - 1; ++y)
        {
          for( size_t z = 1; z < m_Data.GetNumberCols() - 1; ++z)
          {
            double sum_x( 0);
            double sum_y( 0);
            double sum_z( 0);

            for( int i = -1; i < 2; ++i)
            {
              for( int j = -1; j < 2; ++j)
              {
                for( int k = -1; k < 2; ++k)
                {
                  // x gradient
                  sum_x += m_Data( x + i, y + j, z + k) * mask_x( i + 1, j + 1, k + 1);
                  // y gradient
                  sum_y += m_Data( x + i, y + j, z + k) * mask_y( i + 1, j + 1, k + 1);
                  // z gradient
                  sum_z += m_Data( x + i, y + j, z + k) * mask_z( i + 1, j + 1, k + 1);
                }
              }
            }

            edges( x, y, z) = math::Absolute( sum_x) + math::Absolute( sum_y) + math::Absolute( sum_z);
          }
        }
      }

      // update the data
      m_Data = edges;

      // end
      return *this;
    }

    //! @brief balance the intensities
    //! @param LENGTH of cube vertex in which intensities are balanced
    Map &Map::BalanceIntensities( const double LENGTH)
    {
      if( LENGTH <= double( 0))
      {
        return *this;
      }

      const int nr_voxels_x( int( ( LENGTH / 2) / m_CellWidth.X()));
      const int nr_voxels_y( int( ( LENGTH / 2) / m_CellWidth.Y()));
      const int nr_voxels_z( int( ( LENGTH / 2) / m_CellWidth.Z()));

      if( nr_voxels_x == 0 || nr_voxels_y == 0 || nr_voxels_z == 0)
      {
        return *this;
      }
      BCL_MessageStd( "balancing intensities");
      math::Tensor< double> balanced( m_Data);
      const double size( ( 2 * nr_voxels_x + 1) * ( 2 * nr_voxels_y + 1) * ( 2 * nr_voxels_z + 1));
      // iterate over all voxel
      for( size_t x = 0; x < m_Data.GetNumberCols(); ++x)
      {
        for( size_t y = 0; y < m_Data.GetNumberRows(); ++y)
        {
          for( size_t z = 0; z < m_Data.NumberLayers(); ++z)
          {
            double sum( 0);
            for( int i = -nr_voxels_x; i < nr_voxels_x; ++i)
            {
              for( int j = -nr_voxels_y; j < nr_voxels_y; ++j)
              {
                for( int k = -nr_voxels_z; k < nr_voxels_z; ++k)
                {
                  sum += m_Data
                    (
                      ( i + x + m_Data.GetNumberCols()) % m_Data.GetNumberCols(),
                      ( j + y + m_Data.GetNumberRows()) % m_Data.GetNumberRows(),
                      ( k + z + m_Data.NumberLayers()) % m_Data.NumberLayers()
                    );
                }
              }
            }
            balanced( x, y, z) *= size / sum;
          }
        }
      }

      m_Data = balanced;
      BCL_MessageStd( "balancing intensities finished");

      // end
      return *this;
    }

    //! @brief convert Map into Spline
    math::TricubicSpline Map::ConvertDensityToSpline() const
    {
      const math::SplineBorderType border[ 3] = { math::e_FirstDer, math::e_FirstDer, math::e_FirstDer};

      // these vectors are used to input the starting points and the
      // grid width delta of every dimension (x, y, z) into the spline

      //calculate temporary origin
      const linal::Vector3D origin
      (
        m_Origin + linal::Vector3D
        (
          ( m_Index( 0) + 0.5) * m_CellWidth.X(),
          ( m_Index( 1) + 0.5) * m_CellWidth.Y(),
          ( m_Index( 2) + 0.5) * m_CellWidth.Z()
        )
      );
      const double start[] =
      {
        origin.Z(), origin.Y(), origin.X()
      };

      const double delta[] =
      {
        m_CellWidth.Z(), m_CellWidth.Y(), m_CellWidth.X()
      };

      const bool lin_cont[ 3] = { true, true, true};

      //this vector controls the behavior of the spline at the beginning and
      //end of every dimension, only has impact for SplineBorderType FIRSTDER

      //every pair describes the value of the first order derivative at start and end
      const storage::Pair< double, double> first_be[ 3] =
      {
        storage::Pair< double, double>( 1, -1),
        storage::Pair< double, double>( 1, -1),
        storage::Pair< double, double>( 1, -1)
      };

      //convert parameter density map in TriCubicSpline
      math::TricubicSpline densitymap_as_spline;
      densitymap_as_spline.Train( border, start, delta, m_Data, lin_cont, first_be);

      return densitymap_as_spline;
    }

    //! @brief orthogonalize map if angles aren't 90,90,90
    Map Map::OrthogonalizeMap() const
    {
      // if angles are already 90,90,90
      if( m_Angle.Z() == 90 && m_Angle.Y() == 90 && m_Angle.X() == 90)
      {
        // return just the original density map
        return Map
        (
          m_Data, m_Index, m_Intervals, m_Length, m_CellWidth, m_Angle, m_Axis, m_Origin
        );
      }

      // if angles are not already 90,90,90 then orthogonalize the map
      // for this calculate the point that is origin + index of the original map in real space coordinates
      // (this is point in map that is closest to the origin)
      const linal::Vector3D begin_of_map_non_orthog_system
      (
        m_CellWidth.X() * m_Index( 0) + m_Origin.X(),
        m_CellWidth.Y() * m_Index( 1) + m_Origin.Y(),
        m_CellWidth.Z() * m_Index( 2) + m_Origin.Z()
      );

      // begin of the original map in orthogonal system is the same (per definition)
      const linal::Vector3D begin_of_map_orthog_system( begin_of_map_non_orthog_system);

      // calculate the end point of original map (point that is farthest away from begin point and origin)
      const linal::Vector3D end_of_map_non_orthog_system
      (
        m_CellWidth.X() * ( m_Index( 0) + GetDimensions()( 0)) + m_Origin.X(),
        m_CellWidth.Y() * ( m_Index( 1) + GetDimensions()( 1)) + m_Origin.Y(),
        m_CellWidth.Z() * ( m_Index( 2) + GetDimensions()( 2)) + m_Origin.Z()
      );

      // convert the end point coordinates of original map into orthogonal system
      const linal::Vector3D end_of_map_orthog_system
      (
        end_of_map_non_orthog_system.X() + cos( math::Angle::Radian( m_Angle.Y())) * m_CellWidth.Z() * GetDimensions()( 2)
          + cos( math::Angle::Radian( m_Angle.X())) * m_CellWidth.Y() * GetDimensions()( 1),
        end_of_map_non_orthog_system.Y() + cos( math::Angle::Radian( m_Angle.Z())) * m_CellWidth.Z() * GetDimensions()( 2)
          + ( sin( math::Angle::Radian( m_Angle.X())) - 1.0) * m_CellWidth.Y() * GetDimensions()( 1),
        end_of_map_non_orthog_system.Z() + ( sin( math::Angle::Radian( m_Angle.Y())) - 1.0) * m_CellWidth.Z() * GetDimensions()( 2)
          + ( sin( math::Angle::Radian( m_Angle.Z())) - 1.0) * m_CellWidth.Z() * GetDimensions()( 2)
      );

      // calculate the new cell width (make sure it is not larger than the original cell width)
      const linal::Vector3D new_cell_width
      (
        std::abs( sin( math::Angle::Radian( m_Angle.X()))) * m_CellWidth.X(),
        std::abs( sin( math::Angle::Radian( m_Angle.Y()))) * m_CellWidth.Y(),
        std::abs( sin( math::Angle::Radian( m_Angle.Z()))) * m_CellWidth.Z()
      );

      // for debug purposes give out the begin and end point of the original map
      BCL_MessageDbg
      (
        "ATOM   9996  N   ASP B   1     " + util::Format().W( 7).FFP( 3).R()( begin_of_map_orthog_system.X()) +
        " " + util::Format().W( 7).FFP( 3).R()( begin_of_map_orthog_system.Y()) +
        " " + util::Format().W( 7).FFP( 3).R()( begin_of_map_orthog_system.Z()) + "  1.00 38.08           N"
      );
      BCL_MessageDbg
      (
        "ATOM   9996  N   ASP B   1     " + util::Format().W( 7).FFP( 3).R()( end_of_map_orthog_system.X()) +
        " " + util::Format().W( 7).FFP( 3).R()( end_of_map_orthog_system.Y()) +
        " " + util::Format().W( 7).FFP( 3).R()( end_of_map_orthog_system.Z()) + "  1.00 38.08           N"
      );

      // determine the begin and end points of the new map in real space coordinates (these points will be different
      // if any angle is greater 90 degrees)
      double begin_of_new_map_X( begin_of_map_orthog_system.X());
      double begin_of_new_map_Y( begin_of_map_orthog_system.Y());
      double begin_of_new_map_Z( begin_of_map_orthog_system.Z());
      double end_of_new_map_X( end_of_map_orthog_system.X());
      double end_of_new_map_Y( end_of_map_orthog_system.Y());
      double end_of_new_map_Z( end_of_map_orthog_system.Z());

      // if one of the angles is greater than 90, than the new map will be larger in one dimension
      // if angle is greater 90
      if( m_Angle.X() > 90)
      {
        // new begin will be old begin minus length of density map in this direction * cos ( 180 - angle)
        begin_of_new_map_X += cos( math::Angle::Radian( m_Angle.X())) * m_CellWidth.Y() * GetDimensions()( 1);

        // new end will be old end plus length of density map in this direction * cos ( 180 - angle)
        end_of_new_map_X -= cos( math::Angle::Radian( m_Angle.X())) * m_CellWidth.Y() * GetDimensions()( 1);
      }

      // if angle is greater 90
      if( m_Angle.Y() > 90)
      {
        // new begin will be old begin minus length of density map in this direction * cos ( 180 - angle)
        begin_of_new_map_X += cos( math::Angle::Radian( m_Angle.Y())) * m_CellWidth.Z() * GetDimensions()( 2);

        // new end will be old end plus length of density map in this direction * cos ( 180 - angle)
        end_of_new_map_X -= cos( math::Angle::Radian( m_Angle.Y())) * m_CellWidth.Z() * GetDimensions()( 2);
      }

      // if angle is greater 90
      if( m_Angle.Z() > 90)
      {
        // new begin will be old begin minus length of density map in this direction * cos ( 180 - angle)
        begin_of_new_map_Y += cos( math::Angle::Radian( m_Angle.Z())) * m_CellWidth.Z() * GetDimensions()( 2);

        // new end will be old end plus length of density map in this direction * cos ( 180 - angle)
        end_of_new_map_Y -= cos( math::Angle::Radian( m_Angle.Z())) * m_CellWidth.Z() * GetDimensions()( 2);
      }

      // for the new orthogonalized map make sure that it has no origin and all the shift is contained in the index
      // calculate the new index
      int new_index_X( ( begin_of_new_map_X) / new_cell_width.X());
      if( new_index_X > 0)
      {
        new_index_X = floor( new_index_X);
      }
      else if( new_index_X < 0)
      {
        new_index_X = ceil( new_index_X);
      }

      // calculate the new index
      int new_index_Y( ( begin_of_new_map_Y) / new_cell_width.Y());
      if( new_index_Y > 0)
      {
        new_index_Y = floor( new_index_Y);
      }
      else if( new_index_Y < 0)
      {
        new_index_Y = ceil( new_index_Y);
      }

      // calculate the new index
      int new_index_Z( ( begin_of_new_map_Z) / new_cell_width.Z());
      if( new_index_Z > 0)
      {
        new_index_Z = floor( new_index_Z);
      }
      else if( new_index_Z < 0)
      {
        new_index_Z = ceil( new_index_Z);
      }

      // calculate new extension of density map (corresponding to number of cols for the new map)
      size_t new_extension_X( ceil( ( end_of_new_map_X - begin_of_new_map_X) / new_cell_width.X()));

      // calculate new extension of density map (corresponding to number of rows for the new map)
      size_t new_extension_Y( ceil( ( end_of_new_map_Y - begin_of_new_map_Y) / new_cell_width.Y()));

      // calculate new extension of density map (corresponding to number of layers for the new map)
      size_t new_extension_Z( ceil( ( end_of_new_map_Z - begin_of_new_map_Z) / new_cell_width.Z()));

      // set begin of new map
      linal::Vector3D begin_of_new_map
      (
        begin_of_new_map_X,
        begin_of_new_map_Y,
        begin_of_new_map_Z
      );
      // set end of new map
      linal::Vector3D end_of_new_map
      (
        end_of_new_map_X,
        end_of_new_map_Y,
        end_of_new_map_Z
      );

      // for debug purposes give out the begin and end point of the new map
      BCL_MessageDbg
      (
        "ATOM   9997  N   ASP B   1     " + util::Format().W( 7).FFP( 3).R()( begin_of_new_map.X()) +
        " " + util::Format().W( 7).FFP( 3).R()( begin_of_new_map.Y()) +
        " " + util::Format().W( 7).FFP( 3).R()( begin_of_new_map.Z()) + "  1.00 38.08           N"
      );
      BCL_MessageDbg
      (
        "ATOM   9998  N   ASP B   1     " + util::Format().W( 7).FFP( 3).R()( end_of_new_map.X()) +
        " " + util::Format().W( 7).FFP( 3).R()( end_of_new_map.Y()) +
        " " + util::Format().W( 7).FFP( 3).R()( end_of_new_map.Z()) + "  1.00 38.08           N"
      );

      // set new index and extension
      linal::VectorND< int, 3> new_index( new_index_X, new_index_Y, new_index_Z);
      storage::VectorND< 3, size_t> new_extension( new_extension_X, new_extension_Y, new_extension_Z);

      // set the new intervals correctly
      linal::VectorND< int, 3> new_intervals;
      for( size_t i( 0); i < 3; ++i)
      {
        if( new_index( i) <= -int( new_extension( i)))
        {
          new_intervals( i) = math::Absolute( new_index( i));
        }
        else if( new_index( i) >= 0)
        {
          new_intervals( i) = new_index( i) + new_extension( i) - 1;
        }
        else
        {
          new_intervals( i) = new_extension( i) - 1;
        }
      }

      // set the new length
      const linal::Vector3D new_length
      (
        new_cell_width.X() * double( new_intervals( 0)),
        new_cell_width.Y() * double( new_intervals( 1)),
        new_cell_width.Z() * double( new_intervals( 2))
      );

      // set the new origin to 0, 0, 0
      linal::Vector3D new_origin( 0.0, 0.0, 0.0);

      // set the new angles to 90, 90, 90
      linal::Vector3D new_angles( 90, 90, 90);

      // convert the original map into a spline
      math::TricubicSpline old_map_as_spline( ConvertDensityToSpline());

      // create tensor of appropriate size for new data, filled with the minimum value of the original map
      math::Tensor< double> new_data
      (
        new_extension_Z, // 2
        new_extension_Y, // 1
        new_extension_X, // 0
        m_Minimum
      );

      // store dimensions of the original map in a storage::VectorND
      linal::VectorND< int, 3> old_dimensions
      (
        int( m_Data.GetNumberCols()),
        int( m_Data.GetNumberRows()),
        int( m_Data.NumberLayers())
      );

      // iterate over the x direction
      for( size_t itr_X( 0); itr_X < new_extension_X; ++itr_X)
      {
        // calculate real space coordinates of this particular point (for X coordinate)
        double coordinates_orthogonal_system_X( ( new_index_X + int( itr_X) + 0.5) * new_cell_width.X());

        // iterate over the y direction
        for( size_t itr_Y( 0); itr_Y < new_extension_Y; ++itr_Y)
        {
          // calculate real space coordinates of this particular point (for Y coordinate)
          double coordinates_orthogonal_system_Y( ( new_index_Y + int( itr_Y) + 0.5) * new_cell_width.Y());

          // iterate over the z direction
          for( size_t itr_Z( 0); itr_Z < new_extension_Z; ++itr_Z)
          {
            // calculate real space coordinates of this particular point (for Y coordinate)
            double coordinates_orthogonal_system_Z( ( new_index_Z + int( itr_Z) + 0.5) * new_cell_width.Z());

            // calculate real space coordinates of this particular point
            const linal::Vector3D coordinates_orthogonal_system
            (
              coordinates_orthogonal_system_X,
              coordinates_orthogonal_system_Y,
              coordinates_orthogonal_system_Z
            );

            // convert these real space coordinates into the non-orthogonal system
            const linal::Vector3D coordinates_non_orthogonal_system
            (
              coordinates_orthogonal_system.X() - cos( math::Angle::Radian( m_Angle.X())) * m_CellWidth.Y() * itr_Y
                - cos( math::Angle::Radian( m_Angle.Y())) * m_CellWidth.Z() * itr_Z,
              coordinates_orthogonal_system.Y() - cos( math::Angle::Radian( m_Angle.Z())) * m_CellWidth.Z() * itr_Z
                - ( sin( math::Angle::Radian( m_Angle.X())) - 1.0) * m_CellWidth.Y() * itr_Y,
              coordinates_orthogonal_system.Z() - ( sin( math::Angle::Radian( m_Angle.Y())) - 1.0) * m_CellWidth.Z() * itr_Z
                - ( sin( math::Angle::Radian( m_Angle.Z())) - 1.0) * m_CellWidth.Z() * itr_Z
            );

            // calculate the voxel in the old map that this real space coordinate corresponds to
            const linal::VectorND< int, 3> old_map_lookup
            (
              (( coordinates_non_orthogonal_system.X() - m_Origin.X()) / m_CellWidth.X()) - m_Index( 0),
              (( coordinates_non_orthogonal_system.Y() - m_Origin.Y()) / m_CellWidth.Y()) - m_Index( 1),
              (( coordinates_non_orthogonal_system.Z() - m_Origin.Z()) / m_CellWidth.Z()) - m_Index( 2)
            );

            // if you are looking up a coordinate that corresponds to a defined voxel in the original map, look it up
            // if this coordinate is outside the original map, don't do anything, new_data already has the minimum
            // value stored by default
            if
            (
              old_map_lookup( 0) >= 0 && old_map_lookup( 0) < old_dimensions( 0) &&
              old_map_lookup( 1) >= 0 && old_map_lookup( 1) < old_dimensions( 1) &&
              old_map_lookup( 2) >= 0 && old_map_lookup( 2) < old_dimensions( 2)
            )
            {
              // lookup intensity in spline and insert into this density map
              new_data( itr_Z, itr_Y, itr_X) =
                old_map_as_spline.F
                (
                  coordinates_non_orthogonal_system.Z(),
                  coordinates_non_orthogonal_system.Y(),
                  coordinates_non_orthogonal_system.X()
                );

              // this is second possible implementation that is not based on the spline
              // it just looks up the intensity at this particular voxel, no averaging for off-center-voxel value is done
//              new_data( itr_Z, itr_Y, itr_X) = m_Data( old_map_lookup( 2), old_map_lookup( 1), old_map_lookup( 0));
            }
          }
        }
      }

      // return the orthogonalized density map
      return Map
      (
        new_data, new_index, new_intervals, new_length, new_cell_width, new_angles, m_Axis, new_origin
      );
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Map::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "hashmap class that reads in mrc files");
      serializer.AddInitializer
        (
         "index",
         "index of map",
         io::Serialization::GetAgent( &m_Index)
         );
      serializer.AddInitializer
        (
         "intervals",
         "intervals along each axis",
         io::Serialization::GetAgent( &m_Intervals)
         );
      serializer.AddInitializer
        (
         "length",
         "size of unit cell in angstrom",
         io::Serialization::GetAgent( &m_Length)
         );
      serializer.AddInitializer
        (
         "cell width",
         "size of each voxel in angstrom",
         io::Serialization::GetAgent( &m_CellWidth)
         );
      serializer.AddInitializer
        (
         "angle",
         "voxel angle in degrees",
         io::Serialization::GetAgent( &m_Angle)
         );
      serializer.AddInitializer
        (
         "axis",
         "column, row, and section axes",
         io::Serialization::GetAgent( &m_Axis)
         );
      serializer.AddInitializer
        (
         "minimum",
         "minimum intensity",
         io::Serialization::GetAgent( &m_Minimum)
         );
      serializer.AddInitializer
        (
         "maximum",
         "maximum intensity",
         io::Serialization::GetAgent( &m_Maximum)
         );
      serializer.AddInitializer
        (
         "mean",
         "mean intensity",
         io::Serialization::GetAgent( &m_Mean)
         );
      serializer.AddInitializer
        (
         "space group",
         "space group number (0 or 1)",
         io::Serialization::GetAgent( &m_SpaceGroup)
         );
      serializer.AddInitializer
        (
         "number bytes symmetry data",
         "number of bytes used for symmetry data (0 or 80)",
         io::Serialization::GetAgent( &m_NumberBytesSymmetryData)
         );
      // serializer.AddInitializer
      //   (
      //    "extra",
      //    "extra data, 0 by default",
      //    io::Serialization::GetAgent( &m_Extra)
      //    );
      serializer.AddInitializer
        (
         "rmsd",
         "rmsd of density to mean density value",
         io::Serialization::GetAgent( &m_Rmsd)
         );
      serializer.AddInitializer
        (
         "origin",
         "origin of map",
         io::Serialization::GetAgent( &m_Origin)
         );
      serializer.AddInitializer
        (
         "data",
         "density as a tensor that also contains the dimensions of the density map",
         io::Serialization::GetAgent( &m_Data)
         );
      serializer.AddInitializer
        (
         "machine stamp",
         "machine stamp in header",
         io::Serialization::GetAgent( &m_MachineStamp)
         );
      serializer.AddInitializer
        (
         "labels",
         "labels in header",
         io::Serialization::GetAgent( &m_Labels)
         );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator
    //! @param MAP_RHS right hand side density map to be copied
    //! @return Map the changed map
    Map &Map::operator=( const Map &MAP_RHS)
    {
      // same map
      if( &MAP_RHS == this)
      {
        return *this;
      }

      // assign members
      m_Index                   = MAP_RHS.m_Index;
      m_Intervals               = MAP_RHS.m_Intervals;
      m_Length                  = MAP_RHS.m_Length;
      m_CellWidth               = MAP_RHS.m_CellWidth;
      m_Angle                   = MAP_RHS.m_Angle;
      m_Axis                    = MAP_RHS.m_Axis;
      m_Minimum                 = MAP_RHS.m_Minimum;
      m_Maximum                 = MAP_RHS.m_Maximum;
      m_Mean                    = MAP_RHS.m_Mean;
      m_SpaceGroup              = MAP_RHS.m_SpaceGroup;
      m_NumberBytesSymmetryData = MAP_RHS.m_NumberBytesSymmetryData;
      m_Rmsd                    = MAP_RHS.m_Rmsd;
      m_Origin                  = MAP_RHS.m_Origin;
      m_Data                    = MAP_RHS.m_Data;
      m_MachineStamp            = MAP_RHS.m_MachineStamp;
      m_Labels                  = MAP_RHS.m_Labels;

      // copy extra
      std::copy( MAP_RHS.m_Extra, MAP_RHS.m_Extra + s_ExtraSize, m_Extra);

      // end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! Write Header
    std::ostream &Map::WriteHeader( std::ostream &OSTREAM) const
    {
      OSTREAM << "Dimensions in x, y, z:\n" << GetDimensions() << '\n'; // Dimensions of the map in x, y ,z
      OSTREAM << "Index:\n"                 << m_Index      << '\n'; // index of map
      OSTREAM << "Intervals:\n"             << m_Intervals  << '\n'; // intervals along each axis
      OSTREAM << "Length x, y, z:\n"        << m_Length     << '\n'; // size of unitcell in Angstrom
      OSTREAM << "Angle  a, b, c:\n"        << m_Angle      << '\n'; // voxel angle in degrees
      OSTREAM << "Axis( fast - slow):\n"    << m_Axis       << '\n'; // col, row and section axis
      OSTREAM << "Minimum Intensity:\n"     << m_Minimum    << '\n'; // Minimal density value
      OSTREAM << "Maximum Intensity:\n"     << m_Maximum    << '\n'; // Maximal density value
      OSTREAM << "Mean Intensity:\n"        << m_Mean       << '\n'; // Mean density value
      OSTREAM << "SpaceGroup: "             << m_SpaceGroup << '\n';
      OSTREAM << "Number sym data bytes: "  << m_NumberBytesSymmetryData << '\n';
      OSTREAM << "extra data: "             << std::string( m_Extra, s_ExtraSize) << '\n';
      OSTREAM << "RMSD:\n"                  << m_Rmsd         << '\n'; // rms of density to mean density value
      OSTREAM << "Origin x, y, z:\n"        << m_Origin       << '\n'; // origin of map
      OSTREAM << "machine stamp:\n"         << m_MachineStamp << '\n'; // machine stamp
      OSTREAM << "labels\n"                 << m_Labels       << '\n'; // labels

      // end
      return OSTREAM;
    }

    //! read Map from mrc file
    //! http://www2.mrc-lmb.cam.ac.uk/image2000.html
    std::istream &Map::ReadMRC( std::istream &ISTREAM, const size_t EXTENDED_HEADER)
    {
      // get length of file:
      ISTREAM.seekg( 0, std::ios::end);
      const std::istream::pos_type file_length( ISTREAM.tellg());
      std::istream::pos_type file_current;
      ISTREAM.seekg( 0, std::ios::beg);

      // determine if it is necessary to swap bytes when reading chars
      const bool swap( CheckSwap( ISTREAM));

      linal::VectorND< int, 3> dimensions;
      //read dimensions; bit 0-11
      dimensions( 0) = ReadInt( ISTREAM, swap); // NX
      dimensions( 1) = ReadInt( ISTREAM, swap); // NY
      dimensions( 2) = ReadInt( ISTREAM, swap); // NZ

      linal::VectorND< size_t, 3> dimensions_copy( dimensions.Begin(), dimensions.End());

      int mode; // MODE; bit 12-15
      mode = ReadInt( ISTREAM, swap);  //read mode of mrc file ( MODE     data type :
                                       //0 image : signed 8-bit bytes range -128 to 127
                                       //1 image : 16-bit halfwords
                                       //2 image : 32-bit reals
                                       //3 transform : complex 16-bit integers
                                       //4 transform : complex 32-bit reals

      BCL_Assert
      (
        mode == 2,
        "MRC file not written in MODE 2 - BCL is not able to read such a format. Written in MODE " +
        util::Format()( mode) + "\ndimensions read so far:\n" + util::Format()( dimensions) +
        "\nand swap is determined to be: " + util::Format()( swap)
      );

      // read Index; bit 16-27
      m_Index( 0)  = ReadInt( ISTREAM, swap); // NXSTART
      m_Index( 1) = ReadInt( ISTREAM, swap); // NYSTART
      m_Index( 2)  = ReadInt( ISTREAM, swap); // NZSTART

      // read Intervals; bit 28-39
      m_Intervals( 0)  = ReadInt( ISTREAM, swap); // MX
      m_Intervals( 1) = ReadInt( ISTREAM, swap); // MY
      m_Intervals( 2)  = ReadInt( ISTREAM, swap); // MZ

      //read width; bit bit 40-51
      for( double *ptr( m_Length.Begin()), *ptr_end( m_Length.End()); ptr != ptr_end; ++ptr)
      {
        *ptr = ReadFloat( ISTREAM, swap); // CELLA
      }

      //calculate m_CellWidth
      m_CellWidth.X() = m_Length.X() / m_Intervals( 0);
      m_CellWidth.Y() = m_Length.Y() / m_Intervals( 1);
      m_CellWidth.Z() = m_Length.Z() / m_Intervals( 2);

      //read Angle; bit 52 - 63
      for( double *ptr( m_Angle.Begin()), *ptr_end( m_Angle.End()); ptr != ptr_end; ++ptr)
      {
        *ptr = ReadFloat( ISTREAM, swap); // CELLB
      }

      //read Axis; bit 64-75
      m_Axis( 0)  = ReadInt( ISTREAM, swap); // MAPC cols
      m_Axis( 1) = ReadInt( ISTREAM, swap); // MAPC rows
      m_Axis( 2)  = ReadInt( ISTREAM, swap); // MAPC sections

      //read minimal, maximal and mean density value of map / will be recalculated; bit 76-87
      m_Minimum = ReadFloat( ISTREAM, swap); // DMIN
      m_Maximum = ReadFloat( ISTREAM, swap); // DMAX
      m_Mean    = ReadFloat( ISTREAM, swap); // DMEAN

      // additional information
      m_SpaceGroup = ReadInt( ISTREAM, swap); // ISPG; bit 88-91
      m_NumberBytesSymmetryData = ReadInt( ISTREAM, swap); // NSYMBT; bit 92-95
      BCL_Assert( m_NumberBytesSymmetryData == 0 || m_NumberBytesSymmetryData == 80, "number of bytes for symmetry data should either be 0 or 80");
      ISTREAM.read( m_Extra, s_ExtraSize); // bit 96-195

      // read origin of map; bit 196-207
      for( double *ptr( m_Origin.Begin()), *ptr_end( m_Origin.End()); ptr != ptr_end; ++ptr)
      {
        *ptr = ReadFloat( ISTREAM, swap); // ORIGIN
      }

      // read char string MAP; bit 208-211
      char tmp[ 5];
      tmp[ 4] = '\0';
      const char *map_ident( "MAP ");
      for( int i( 0); i < 4; ++i)
      {
        *( tmp + i) = ISTREAM.get();
      }
      if( std::string( tmp).compare( map_ident) != 0)
      {
        BCL_MessageCrt( "given mrc map does not contain \"MAP\" at pos 53 in header but: " + std::string( tmp));
      }

      // read machine stamp; bit 212-215
      m_MachineStamp = ReadInt( ISTREAM, swap);

      //checks if m_Axis is a right handed coordinate system - and orders if necessary
      if( m_Axis( 0) != 1 || m_Axis( 1) != 2)
      {
        BCL_MessageCrt
        (
          "Axis have a different order than x, y, z - all information will be ordered x, y, z. No Guarantee that this works!"
        );

        // apparently intervals, length , cellwidth and angle don't have to be ordered along with the axes
        // the index definitely has to be ordered
        // origin probably has to be ordered, but example that we used contained origin (0, 0, 0), so we couldn't
        // check whether origin has to be ordered...
        dimensions  = OrderValues( m_Axis, dimensions);
        m_Index     = OrderValues( m_Axis, m_Index);
//        m_Intervals = OrderValues( m_Axis, m_Intervals); // already ordered
//        m_Length    = OrderValues( m_Axis, m_Length);    // already ordered
//        m_CellWidth = OrderValues( m_Axis, m_CellWidth); // already ordered
//        m_Angle     = OrderValues( m_Axis, m_Angle);     // already ordered
        m_Origin    = OrderValues( m_Axis, m_Origin);
      }

      m_Data = math::Tensor< double>( dimensions( 2), dimensions( 1), dimensions( 0));

      //read rmsd; bit 216-219
      m_Rmsd = ReadFloat( ISTREAM, swap);

      // read number labels
      int number_labels;
      number_labels = ReadInt( ISTREAM, swap); // NLABL; bit 220-223
      BCL_Assert( number_labels <= s_MaxNumberLabels, "there should not be more than 10 labels, but there are: " + util::Format()( number_labels));

      // remove all current labels
      m_Labels.Reset();
      // read all labels
      for( int i( 0); i < number_labels; ++i)
      {
        char label[ s_LabelLength];
        ISTREAM.read( label, s_LabelLength);
        m_Labels.PushBack( std::string( label, s_LabelLength));
      }

      // read density into the tensor
      ISTREAM.seekg( 4 * 256 + EXTENDED_HEADER); //jump to end of header
      file_current = ISTREAM.tellg();

      if( int( file_length - file_current) != int( 4 * m_Data.GetSize()))
      {
        const std::istream::pos_type additional_header( ( file_length - file_current) - int( 4 * m_Data.GetSize()));
        std::ostringstream msg;
        msg << "header does not have a size of 1024 bits!\n"
            << "size of extended header is: " << additional_header << '\n'
            << "proceeding at logical position according to header information for number of stored intensities = "
            << GetSize() << '\n';
        BCL_MessageVrb( msg.str());
        ISTREAM.seekg( file_current + additional_header);
      }

      //read all intensity values
      if( m_Axis( 0) != 1 || m_Axis( 1) != 2)
      {
        linal::VectorND< size_t, 3> index;
        for( index( 2) = 0; index( 2) < dimensions_copy( 2); ++index( 2))
        {
          for( index( 1) = 0; index( 1) < dimensions_copy( 1); ++index( 1))
          {
            for( index( 0) = 0; index( 0) < dimensions_copy( 0); ++index( 0))
            {
              operator()( OrderValues( m_Axis, index)) = ReadFloat( ISTREAM, swap);
            }
          }
        }

        m_Axis = OrderValues( m_Axis, m_Axis);
      }
      else
      {
        for( double *ptr( m_Data.Begin()), *ptr_end( m_Data.End()); ptr != ptr_end && !( ISTREAM.eof()); ++ptr)
        {
          *ptr = ReadFloat( ISTREAM, swap);
        }
      }

      // calculate the min max mean and rmsd
      CalculateMinMaxMeanRmsd();

      // only use index if origin was given, but index is 0
      if( m_Index( 0) == 0 && m_Index( 1) == 0 && m_Index( 2) == 0)
      {
        m_Index( 0) = int( m_Origin( 0) / m_CellWidth( 0));
        m_Index( 1) = int( m_Origin( 1) / m_CellWidth( 1));
        m_Index( 2) = int( m_Origin( 2) / m_CellWidth( 2));
        m_Origin = 0.0;
      }
      else // check if index is given, that origin is not used
      {
        if( !math::EqualWithinTolerance( 0.0, m_Origin.SquareNorm()))
        {
          BCL_MessageCrt
          (
            "map read that contains an origin and an index:\n" + util::Format()( m_Origin) + "\n" +
            util::Format()( m_Index) + "\norigin will be set to zero, assuming that the origin = index * voxel size"
          );
          m_Origin = 0.0;
        }
      }

      // end
      return ISTREAM;
    }

    //! write Map to mrc file
    //! http://www2.mrc-lmb.cam.ac.uk/image2000.html
    std::ostream &Map::WriteMRC( std::ostream &OSTREAM) const
    {
      int tmpint;
      float tmpfloat;
      tmpint = m_Data.GetNumberCols();       //write length in x direction;
      OSTREAM.write( reinterpret_cast< const char *>( &tmpint), 4);
      tmpint = m_Data.GetNumberRows();       //write length in y direction
      OSTREAM.write( reinterpret_cast< const char *>( &tmpint), 4);
      tmpint = m_Data.NumberLayers();     //write length in z direction;
      OSTREAM.write( reinterpret_cast< const char *>( &tmpint), 4);

      tmpint = 2;
      //write mode of mrc file MODE     data type :
      //0 image : signed 8-bit bytes range -128 to 127
      //1 image : 16-bit halfwords
      //2 image : 32-bit reals
      //3 transform : complex 16-bit integers
      //4 transform : complex 32-bit reals
      OSTREAM.write( reinterpret_cast< const char *>( &tmpint), 4);

      //write Index
      OSTREAM.write( reinterpret_cast< const char *>( &m_Index( 0)), 4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_Index( 1)), 4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_Index( 2)), 4);

      //write Intervals
      OSTREAM.write( reinterpret_cast< const char *>( &m_Intervals( 0)), 4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_Intervals( 1)), 4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_Intervals( 2)), 4);

      //write width
      for( const double *ptr( m_Length.Begin()), *ptr_end( m_Length.End()); ptr != ptr_end; ++ptr)
      {
        tmpfloat = float( *ptr);
        OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);
      }

      //write Angle
      for( const double *ptr( m_Angle.Begin()), *ptr_end( m_Angle.End()); ptr != ptr_end; ++ptr)
      {
        tmpfloat = float( *ptr);
        OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);
      }

      //write Axis
      OSTREAM.write( reinterpret_cast< const char *>( &m_Axis( 0)),  4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_Axis( 1)), 4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_Axis( 2)),  4);

      //write minimal, maximal and mean density value of map
      tmpfloat = float( m_Minimum);
      OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);
      tmpfloat = float( m_Maximum);
      OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);
      tmpfloat = float( m_Mean);
      OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);

      // additional information
      OSTREAM.write( reinterpret_cast< const char *>( &m_SpaceGroup), 4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_NumberBytesSymmetryData), 4);
      OSTREAM.write( m_Extra, s_ExtraSize);

      //write origin of map
      for( const double *ptr( m_Origin.Begin()), *ptr_end( m_Origin.End()); ptr != ptr_end; ++ptr)
      {
        tmpfloat = float( *ptr);
        OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);
      }

      // insert MAP
      const char *map_ident( "MAP ");
      OSTREAM.write( map_ident, 4);

      // insert machine stamp
      OSTREAM.write( reinterpret_cast< const char *>( &m_MachineStamp), 4);

      // write rms
      tmpfloat = float( m_Rmsd);
      OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);

      // write number labels
      const int number_labels( m_Labels.GetSize());
      OSTREAM.write( reinterpret_cast< const char *>( &number_labels), 4);

      // write all labels
      for( int i( 0); i < number_labels; ++i)
      {
        BCL_Assert( int( m_Labels( i).size()) == s_LabelLength, "label dose not have the required length!");
        OSTREAM.write( m_Labels( i).c_str(), s_LabelLength);
      }

      // fill stream with empty labels
      OSTREAM.write( std::string( ( s_MaxNumberLabels - number_labels) * s_LabelLength, ' ').c_str(), ( s_MaxNumberLabels - number_labels) * s_LabelLength);

      // write density from the tensor
      for( const double *ptr( m_Data.Begin()), *ptr_end( m_Data.End()); ptr != ptr_end; ++ptr)
      {
        tmpfloat = float( *ptr);
        OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);
      }

      // end
      return OSTREAM;
    }

    //! write Map to std::ostream using the given util::Format
    std::ostream &Map::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Index,     OSTREAM, INDENT) << '\n'; // index of map
      io::Serialize::Write( m_Intervals, OSTREAM, INDENT) << '\n'; // intervals along each axis
      io::Serialize::Write( m_Length,    OSTREAM, INDENT) << '\n'; // size of unitcell in Angstrom
      io::Serialize::Write( m_CellWidth, OSTREAM, INDENT) << '\n'; // size of each voxel in Angstroem
      io::Serialize::Write( m_Angle,     OSTREAM, INDENT) << '\n'; // voxel angle in degrees
      io::Serialize::Write( m_Axis,      OSTREAM, INDENT) << '\n'; // col, row and section axis
      io::Serialize::Write( m_Minimum,   OSTREAM, INDENT) << '\n'; // Minimal density value
      io::Serialize::Write( m_Maximum,   OSTREAM, INDENT) << '\n'; // Maximal density value
      io::Serialize::Write( m_Mean,      OSTREAM, INDENT) << '\n'; // Mean density value
      io::Serialize::Write( m_Rmsd,      OSTREAM, INDENT) << '\n'; // rmsd of density to mean density value
      io::Serialize::Write( m_Origin,    OSTREAM, INDENT) << '\n'; // origin of map

      // density as a tensor does also contain the dimensions of the density map
      io::Serialize::Write( m_Data,      OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! read Map from io::IFStream
    std::istream &Map::Read( std::istream &ISTREAM)
    {
      //read members
      io::Serialize::Read( m_Index,     ISTREAM); // index of map
      io::Serialize::Read( m_Intervals, ISTREAM); // intervals along each axis
      io::Serialize::Read( m_Length,    ISTREAM); // size of unitcell in Angstrom
      io::Serialize::Read( m_CellWidth, ISTREAM); // size of each voxel in Angstroem
      io::Serialize::Read( m_Angle,     ISTREAM); // voxel angle in degrees
      io::Serialize::Read( m_Axis,      ISTREAM); // col, row and section axis
      io::Serialize::Read( m_Minimum,   ISTREAM); // Minimal density value
      io::Serialize::Read( m_Maximum,   ISTREAM); // Maximal density value
      io::Serialize::Read( m_Mean,      ISTREAM); // Mean density value
      io::Serialize::Read( m_Rmsd,      ISTREAM); // rmsd of density to mean density value
      io::Serialize::Read( m_Origin,    ISTREAM); // origin of map
      io::Serialize::Read( m_Data,      ISTREAM); // density as a tensor does also contain the dimensions of the density map

      // end
      return ISTREAM;
    }

    //!checks stream of mrcfile wether the bytes are swapped
    bool Map::CheckSwap( std::istream &ISTREAM)
    {
      unsigned char *cptr;
      int unswapped, swapped;

      //read first number in header
      ISTREAM.read( reinterpret_cast< char *>( &unswapped), 4);
      //set pointer on stream back to begin of stream
      ISTREAM.seekg( 0);
      swapped = unswapped;

      //swap byte 0, 3 and 1, 2
      cptr = ( unsigned char *)&swapped;
      std::swap( cptr[0], cptr[3]);
      std::swap( cptr[1], cptr[2]);

      BCL_MessageDbg
      (
        "unswapped char: " + util::Format()( unswapped) + " swapped char: " + util::Format()( swapped)
      );

      //if swapped is smaller unswapped then all bytes of the file have to be swapped
      return std::abs( ( long int)( swapped)) < std::abs( ( long int)( unswapped));
    }

    //!read int from istream and swap if necessary
    int Map::ReadInt( std::istream &ISTREAM, const bool SWAP)
    {
      unsigned char *cptr;
      int tmp;

      ISTREAM.read( reinterpret_cast< char *>( &tmp), 4);
      if( SWAP)
      {
        cptr = ( unsigned char *)&tmp;
        std::swap( cptr[0], cptr[3]);
        std::swap( cptr[1], cptr[2]);
      }
      return tmp;
    }

    //!read float from istream and swap if necessary
    float Map::ReadFloat( std::istream &ISTREAM, const bool SWAP)
    {
      unsigned char *cptr;
      float tmp;

      ISTREAM.read( reinterpret_cast< char *>( &tmp), 4);
      if( SWAP)
      {
        cptr = ( unsigned char *)&tmp;
        std::swap( cptr[0], cptr[3]);
        std::swap( cptr[1], cptr[2]);
      }
      return tmp;
    }

    //!order values according to given Axis order
    linal::Vector3D Map::OrderValues( const linal::VectorND< int, 3> &AXIS, const linal::Vector3D &VALUES)
    {
      storage::VectorND< 3, size_t> index_xyz;
      index_xyz( AXIS( 0) - 1) = 0;
      index_xyz( AXIS( 1) - 1) = 1;
      index_xyz( AXIS( 2) - 1) = 2;

      linal::Vector3D new_values( VALUES);
      new_values( 0) = VALUES( index_xyz( 0));
      new_values( 1) = VALUES( index_xyz( 1));
      new_values( 2) = VALUES( index_xyz( 2));

      return new_values;
    }

  } // namespace density
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
#include "density/bcl_density_map_cylindrical.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_cylinder_coordinates.h"
#include "density/bcl_density_map.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_histogram_2d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MapCylindrical::s_Instance
    (
      GetObjectInstances().AddInstance( new MapCylindrical())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    MapCylindrical::MapCylindrical() :
      m_Body(),
      m_HeightResolution(),
      m_RadiusResolution(),
      m_AngleResolution(),
      m_Data(),
      m_Minimum(),
      m_Maximum(),
      m_Mean()
    {
    }

    //! construct MapCylindrical from given Parameters
    //! @param BODY main axis of cylindrical density map as linesegment3D
    //! @param DENSITY_MAP the original euclidean density map from which the cylindrical density map is constructed
    //! @param HEIGHT_RESOLUTION "voxel size" in height direction
    //! @param RADIUS_RESOLUTION "voxel size" in radius direction
    //! @param NUMBER_WEDGES number of wedges that the density is divided into
    //! @param UPPER_RADIUS the maximal radius around main axis that cylindrical density map extents to
    MapCylindrical::MapCylindrical
    (
      const assemble::SSEGeometryInterface &BODY,
      const Map &DENSITY_MAP,
      const double &HEIGHT_RESOLUTION,
      const double &RADIUS_RESOLUTION,
      const size_t &NUMBER_WEDGES,
      const double &UPPER_RADIUS
    ) :
      m_Body( BODY),
      m_HeightResolution( HEIGHT_RESOLUTION),
      m_RadiusResolution( RADIUS_RESOLUTION),
      m_AngleResolution( 2 * math::g_Pi / NUMBER_WEDGES),
      m_Data(),
      m_Minimum(),
      m_Maximum(),
      m_Mean()
    {
      BCL_Assert( BODY.GetOrientation().IsDefined(), "Body has to be defined");
      BCL_Assert
      (
        HEIGHT_RESOLUTION > double( 0) && RADIUS_RESOLUTION > double( 0) && NUMBER_WEDGES > 0,
        "height and radius resolution have to be larger than 0, number of wedges has to be larger than 0"
      );
      BCL_Assert( UPPER_RADIUS > RADIUS_RESOLUTION, "maximal radius has to be greater than radius resolution");

      //convert parameter density map in TriCubicSpline
      math::TricubicSpline densitymap_as_spline( DENSITY_MAP.ConvertDensityToSpline());

      //calculate the number of voxel in height, radius and angle
      const size_t number_height_voxels( size_t( 2 * BODY.GetExtent( coord::GetAxes().e_Z) / HEIGHT_RESOLUTION));
      const size_t number_radius_voxels( size_t( UPPER_RADIUS / RADIUS_RESOLUTION));
      const size_t number_angle_voxels( NUMBER_WEDGES);

      // resize data (create tensor of appropriate size, filled with zeros)
      m_Data =
        math::Tensor< double>
       (
         number_height_voxels,
         number_radius_voxels,
         number_angle_voxels,
         double( 0)
       );

      // iterate over height in resolution steps
      for( size_t itr_length( 0); itr_length < number_height_voxels; ++itr_length)
      {
        // iterate over radius in resolution steps
        for( size_t itr_radius( 0); itr_radius < number_radius_voxels; ++itr_radius)
        {
          // iterate over angles in resolution steps
          for( size_t itr_angle( 0); itr_angle < number_angle_voxels; ++itr_angle)
          {
            // create cylindrical coordinates
            const coord::CylinderCoordinates cylind_coord
            (
              ( int( itr_length) - int( number_height_voxels / 2)) * m_HeightResolution,
              itr_radius * m_RadiusResolution,
              itr_angle * m_AngleResolution
            );

            // convert to Cartesian coordinates
            linal::Vector3D cart_coord( cylind_coord.GetCartesianCoordinates());

            // transform Cartesian coordinates with Body
            cart_coord.Transform( BODY.GetOrientation());

            // lookup intensity in spline and insert into this density map
            m_Data( itr_length, itr_radius, itr_angle) =
              densitymap_as_spline.F( cart_coord.Z(), cart_coord.Y(), cart_coord.X());
          }
        }
      }
    }

    //! construct MapCylindrical from given Parameters
    //! @param BODY
    //! @param SPLINE the spline of the original euclidean density map
    //! @param HEIGHT_RESOLUTION "voxel size" in height direction
    //! @param RADIUS_RESOLUTION "voxel size" in radius direction
    //! @param NUMBER_WEDGES number of wedges that the density is divided into
    //! @param UPPER_RADIUS the maximal radius around main axis that cylindrical density map extents to
    MapCylindrical::MapCylindrical
    (
      const assemble::SSEGeometryInterface &BODY,
      math::TricubicSpline &SPLINE,
      const double &HEIGHT_RESOLUTION,
      const double &RADIUS_RESOLUTION,
      const size_t &NUMBER_WEDGES,
      const double &UPPER_RADIUS
    ) :
      m_Body( BODY),
      m_HeightResolution( HEIGHT_RESOLUTION),
      m_RadiusResolution( RADIUS_RESOLUTION),
      m_AngleResolution( 2 * math::g_Pi / NUMBER_WEDGES),
      m_Data(),
      m_Minimum(),
      m_Maximum(),
      m_Mean()
    {
      BCL_Assert( BODY.GetOrientation().IsDefined(), "Body has to be defined");
      BCL_Assert
      (
        HEIGHT_RESOLUTION > double( 0) && RADIUS_RESOLUTION > double( 0) && NUMBER_WEDGES > 0,
        "height and radius resolution have to be larger than 0, number of wedges has to be larger than 0"
      );
      BCL_Assert( UPPER_RADIUS > RADIUS_RESOLUTION, "maximal radius has to be greater than radius resolution");

      //calculate the number of voxel in height, radius and angle
      const size_t number_height_voxels( size_t( 2 * BODY.GetExtent( coord::GetAxes().e_Z) / HEIGHT_RESOLUTION));
      const size_t number_radius_voxels( size_t( UPPER_RADIUS / RADIUS_RESOLUTION));
      const size_t number_angle_voxels( NUMBER_WEDGES);

      // resize data (create tensor of appropriate size, filled with zeros)
      m_Data =
        math::Tensor< double>
       (
         number_height_voxels,
         number_radius_voxels,
         number_angle_voxels,
         double( 0)
       );

      // iterate over height in resolution steps
      for( size_t itr_length( 0); itr_length < number_height_voxels; ++itr_length)
      {
        // iterate over radius in resolution steps
        for( size_t itr_radius( 0); itr_radius < number_radius_voxels; ++itr_radius)
        {
          // iterate over angles in resolution steps
          for( size_t itr_angle( 0); itr_angle < number_angle_voxels; ++itr_angle)
          {
            // create cylindrical coordinates
            const coord::CylinderCoordinates cylind_coord
            (
              ( int( itr_length) - int( number_height_voxels / 2)) * m_HeightResolution,
              itr_radius * m_RadiusResolution,
              itr_angle * m_AngleResolution
            );

            // convert to Cartesian coordinates
            linal::Vector3D cart_coord( cylind_coord.GetCartesianCoordinates());

            // transform Cartesian coordinates with Body
            cart_coord.Transform( BODY.GetOrientation());

            // lookup intensity in spline and insert into this density map
            m_Data( itr_length, itr_radius, itr_angle) =
              SPLINE.F( cart_coord.Z(), cart_coord.Y(), cart_coord.X());
          }
        }
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate list of cylindrical density maps from list of bodies (helices) and density map
    //! @param BODIES list of bodies (helices)
    //! @param DENSITY_MAP Cartesian density map
    //! @param HEIGHT_RESOLUTION "voxel size" in height direction
    //! @param RADIUS_RESOLUTION "voxel size" in radius direction
    //! @param NUMBER_WEDGES number of wedges that the density is divided into
    //! @param UPPER_RADIUS the maximal radius around main axis that cylindrical density map extents to
    //! @return return a list of cylindrical density maps
    storage::List< MapCylindrical> MapCylindrical::CalculateCylindricalMaps
    (
      const util::SiPtrList< const assemble::SSEGeometryInterface> &BODIES,
      const Map &DENSITY_MAP,
      const double &HEIGHT_RESOLUTION,
      const double &RADIUS_RESOLUTION,
      const size_t &NUMBER_WEDGES,
      const double &UPPER_RADIUS
    )
    {
      // convert density map into spline
      math::TricubicSpline densitymap_as_spline( DENSITY_MAP.ConvertDensityToSpline());

      // Initialize empty list of cylindrical density maps
      storage::List< MapCylindrical> list_of_cylindrical_density_maps;

      // iterate over all bodies
      for
      (
        util::SiPtrList< const assemble::SSEGeometryInterface>::const_iterator list_itr( BODIES.Begin()),
          list_itr_end( BODIES.End());
        list_itr != list_itr_end;
        ++list_itr
      )
      {
        // use body and spline to construct cylindrical density map
        // insert cylindrical density map into our list of cylindrical density maps
        list_of_cylindrical_density_maps.PushBack
        (
          MapCylindrical
          (
            **list_itr,
            densitymap_as_spline,
            HEIGHT_RESOLUTION,
            RADIUS_RESOLUTION,
            NUMBER_WEDGES,
            UPPER_RADIUS
          )
        );
      }

      // return list of cylindrical density maps
      return list_of_cylindrical_density_maps;
    }

    //! @brief add up all wedges to get a 2D profile (height vs radius - with intensity coloring)
    //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
    //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
    //! @return 2D intensity profile along main axis (in histogram x is height, y is radius)
    math::Histogram2D MapCylindrical::TwoDProfileHeightRadius( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const
    {

      BCL_Assert( LOWER_RADIUS > double( 0) && LOWER_RADIUS < UPPER_RADIUS, "minimal radius must be between 0 and upper radius");

      // initialize Histogram2D of appropriate size to store intensities for all height/radius pairs
      math::Histogram2D two_d_profile
      (
        storage::VectorND< 2, double>( 0.0, LOWER_RADIUS), // min height and radius
        storage::VectorND< 2, double>( m_HeightResolution, m_RadiusResolution), // bin size for height and radius
        storage::VectorND< 2, size_t>( m_Data.NumberLayers(), ( ( UPPER_RADIUS - LOWER_RADIUS) / m_RadiusResolution)) // number bins for height and radius
      );
      // initialize minimal and maximal radii for iteration (in terms of Angstroem and bins)
      const size_t lower_radius_itr( size_t( LOWER_RADIUS / m_RadiusResolution));
      const size_t upper_radius_itr( size_t( UPPER_RADIUS / m_RadiusResolution));

      // iterate over height in resolution steps
      for( size_t itr_length( 0); itr_length < m_Data.NumberLayers(); ++itr_length)
      {
        // iterate over radius in resolution steps
        for( size_t itr_radius( lower_radius_itr); itr_radius < upper_radius_itr; ++itr_radius)
        {
          // initialize temporary sum of intensities over all angles to 0.0
          double sum_intensities( 0.0);

          // iterate over angles in resolution steps
          // optimize this by using statistics::Sum
          for( size_t itr_angle( 0); itr_angle < m_Data.GetNumberCols(); ++itr_angle)
          {
            sum_intensities += m_Data( itr_length, itr_radius, itr_angle);
          }

          // initialize vector of correct height/radius coordinates to pushback into 2D profile
          const storage::VectorND< 2, double> values_for_profile
          (
            ( itr_length + 0.5) * m_HeightResolution,
            ( itr_radius + 0.5) * m_RadiusResolution
          );
          // pushback appropriate sum of intensities over the angles into 2D profile
          two_d_profile.PushBack( values_for_profile, sum_intensities);
        }
      }

      //return 2D profile
      return two_d_profile;
    }

    //! @brief add up all radii to get a 2D profile (height vs angle - with intensity coloring)
    //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
    //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
    //! @return 2D intensity profile along main axis (in histogram x is height, y is angle)
    math::Histogram2D MapCylindrical::TwoDProfileHeightAngle( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const
    {

      BCL_Assert( LOWER_RADIUS > double( 0) && LOWER_RADIUS < UPPER_RADIUS, "minimal radius must be between 0 and upper radius");

      // initialize Histogram2D of appropriate size to store intensities for all height/angle pairs
      math::Histogram2D two_d_profile
      (
        storage::VectorND< 2, double>( 0.0, 0.0), // min height and angle
        storage::VectorND< 2, double>( m_HeightResolution, m_AngleResolution), // binsize for height and angle
        storage::VectorND< 2, size_t>( m_Data.NumberLayers(), m_Data.GetNumberCols()) // number bins for height and angle
      );

      // initialize minimal and maximal radii for iteration (in terms of Angstroem and bins)
      const size_t lower_radius_itr( size_t( LOWER_RADIUS / m_RadiusResolution));
      const size_t upper_radius_itr( size_t( UPPER_RADIUS / m_RadiusResolution));

      // iterate over height in resolution steps
      for( size_t itr_length( 0); itr_length < m_Data.NumberLayers(); ++itr_length)
      {
        // iterate over angles in resolution steps
        for( size_t itr_angle( 0); itr_angle < m_Data.GetNumberCols(); ++itr_angle)
        {
          // initialize temporary sum of intensities over all radii to 0.0
          double sum_intensities( 0.0);

          // iterate over radius in resolution steps (exclude innermost voxels that belong to density rod itself)
          // optimize this by using statistics::Sum
          for
          (
            size_t itr_radius( lower_radius_itr);
            itr_radius < upper_radius_itr;
            ++itr_radius
          )
          {
            sum_intensities += m_Data( itr_length, itr_radius, itr_angle);
          }

          // initialize vector of correct height/angle coordinates to pushback into 2D profile
          const storage::VectorND< 2, double> values_for_profile
          (
            ( itr_length + 0.5) * m_HeightResolution,
            ( itr_angle + 0.5) * m_AngleResolution
          );
          // pushback appropriate sum of intensities over the angles into 2D profile
          two_d_profile.PushBack( values_for_profile, sum_intensities);
        }
      }

      //return 2D profile
      return two_d_profile;
    }

    //! @brief add up all layers to get a 2D profile (radius vs angle - with intensity coloring)
    //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
    //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
    //! @return 2D intensity profile along main axis (in histogram radius is height, y is angle)
    math::Histogram2D MapCylindrical::TwoDProfileRadiusAngle( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const
    {

      BCL_Assert( LOWER_RADIUS > double( 0) && LOWER_RADIUS < UPPER_RADIUS, "minimal radius must be between 0 and upper radius");

      // initialize Histogram2D of appropriate size to store intensities for all radius/angle pairs
      math::Histogram2D two_d_profile
      (
        storage::VectorND< 2, double>( LOWER_RADIUS, 0.0), // min radius and angle
        storage::VectorND< 2, double>( m_RadiusResolution, m_AngleResolution), // bin size for radius and angle
        storage::VectorND< 2, size_t>( ( ( UPPER_RADIUS - LOWER_RADIUS) / m_RadiusResolution), m_Data.GetNumberCols()) // number bins for radius and angle
      );
      // initialize minimal and maximal radii for iteration (in terms of Angstroem and bins)
      const size_t lower_radius_itr( size_t( LOWER_RADIUS / m_RadiusResolution));
      const size_t upper_radius_itr( size_t( UPPER_RADIUS / m_RadiusResolution));

      // iterate over radius in resolution steps
      for( size_t itr_radius( lower_radius_itr); itr_radius < upper_radius_itr; ++itr_radius)
      {
        // iterate over angle in resolution steps
        for( size_t itr_angle( 0); itr_angle < m_Data.GetNumberCols(); ++itr_angle)
        {
          // initialize temporary sum of intensities over all angles to 0.0
          double sum_intensities( 0.0);

          // iterate over height in resolution steps
          // optimize this by using statistics::Sum
          for( size_t itr_length( 0); itr_length < m_Data.NumberLayers(); ++itr_length)
          {
            sum_intensities += m_Data( itr_length, itr_radius, itr_angle);
          }

          // initialize vector of correct radius/angle coordinates to pushback into 2D profile
          const storage::VectorND< 2, double> values_for_profile
          (
            ( itr_radius + 0.5) * m_RadiusResolution,
            ( itr_angle + 0.5) * m_AngleResolution
          );
          // pushback appropriate sum of intensities over the height into 2D profile
          two_d_profile.PushBack( values_for_profile, sum_intensities);
        }
      }

      //return 2D profile
      return two_d_profile;
    }

    //! @brief add up all wedges and radii to get a 1D profile (height vs intensity)
    //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
    //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
    //! @return 1D intensity profile along main axis
    math::Histogram MapCylindrical::OneDProfileHeight( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const
    {

      BCL_Assert( LOWER_RADIUS > double( 0) && LOWER_RADIUS < UPPER_RADIUS, "minimal radius must be between 0 and upper radius");

      // initialize Histogram of appropriate size to store intensities for all height values
      const double min_height( 0.0);
      const double binsize_height( m_HeightResolution);
      const size_t number_bins_height( m_Data.NumberLayers());
      math::Histogram one_d_profile( min_height, binsize_height, number_bins_height);

      // initialize minimal and maximal radii for iteration (in terms of Angstroem and bins)
      const size_t lower_radius_itr( size_t( LOWER_RADIUS / m_RadiusResolution));
      const size_t upper_radius_itr( size_t( UPPER_RADIUS / m_RadiusResolution));

      // iterate over height in resolution steps
      for( size_t itr_length( 0); itr_length < m_Data.NumberLayers(); ++itr_length)
      {
        // initialize temporary sum of intensities over all radii and angles to 0.0
        double sum_intensities( 0.0);

        // iterate over radius in resolution steps (exclude innermost voxels that belong to density rod itself)
        for( size_t itr_radius( lower_radius_itr); itr_radius < upper_radius_itr; ++itr_radius)
        {
          // iterate over angles in resolution steps
          for( size_t itr_angle( 0); itr_angle < m_Data.GetNumberCols(); ++itr_angle)
          {
            sum_intensities += m_Data( itr_length, itr_radius, itr_angle);
          }

          // initialize vector of correct height coordinates to pushback into 1D profile
          const double value_for_profile
          (
            ( itr_length + 0.5) * m_HeightResolution
          );
          // pushback appropriate sum of intensities over the radii and angles into 1D profile
          one_d_profile.PushBack( value_for_profile, sum_intensities);
        }
      }

      //return 1D profile
      return one_d_profile;
    }

    //! @brief add up all wedges and layers to get a 1D profile (radius vs intensity)
    //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
    //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
    //! @return 1D intensity profile along main axis
    math::Histogram MapCylindrical::OneDProfileRadius( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const
    {

      BCL_Assert( LOWER_RADIUS > double( 0) && LOWER_RADIUS < UPPER_RADIUS, "minimal radius must be between 0 and upper radius");

      // initialize Histogram of appropriate size to store intensities for all radius values
      const double min_radius( LOWER_RADIUS);
      const double binsize_radius( m_RadiusResolution);
      const size_t number_bins_radius( ( UPPER_RADIUS - LOWER_RADIUS) / m_RadiusResolution);
      math::Histogram one_d_profile( min_radius, binsize_radius, number_bins_radius);

      // initialize minimal and maximal radii for iteration (in terms of Angstroem and bins)
      const size_t lower_radius_itr( size_t( LOWER_RADIUS / m_RadiusResolution));
      const size_t upper_radius_itr( size_t( UPPER_RADIUS / m_RadiusResolution));

      // iterate over radius in resolution steps
      for( size_t itr_radius( lower_radius_itr); itr_radius < upper_radius_itr; ++itr_radius)
      {
        // initialize temporary sum of intensities over all layers and angles to 0.0
        double sum_intensities( 0.0);

        // iterate over height in resolution steps
        for( size_t itr_length( 0); itr_length < m_Data.NumberLayers(); ++itr_length)
        {
          // iterate over angles in resolution steps
          for( size_t itr_angle( 0); itr_angle < m_Data.GetNumberCols(); ++itr_angle)
          {
            sum_intensities += m_Data( itr_length, itr_radius, itr_angle);
          }

          // initialize vector of correct radius coordinates to pushback into 1D profile
          const double value_for_profile
          (
            ( itr_radius + 0.5) * m_RadiusResolution
          );
          // pushback appropriate sum of intensities over the height and angles into 1D profile
          one_d_profile.PushBack( value_for_profile, sum_intensities);
        }
      }

      //return 1D profile
      return one_d_profile;
    }

    //! @brief add up all wedges and layers to get a 1D profile (angle vs intensity)
    //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
    //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
    //! @return 1D intensity profile along main axis
    math::Histogram MapCylindrical::OneDProfileAngle( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const
    {
      // initialize Histogram of appropriate size to store intensities for all angle values
      const double min_angle( 0.0);
      const double binsize_angle( m_AngleResolution);
      const size_t number_bins_angle( m_Data.GetNumberCols());
      math::Histogram one_d_profile( min_angle, binsize_angle, number_bins_angle);

      // initialize minimal and maximal radii for iteration (in terms of Angstroem and bins)
      const size_t lower_radius_itr( size_t( LOWER_RADIUS / m_RadiusResolution));
      const size_t upper_radius_itr( size_t( UPPER_RADIUS / m_RadiusResolution));

      // iterate over angle in resolution steps
      for( size_t itr_angle( 0); itr_angle < m_Data.GetNumberCols(); ++itr_angle)
      {
        // initialize temporary sum of intensities over all radii and height to 0.0
        double sum_intensities( 0.0);

        // iterate over radius in resolution steps (exclude innermost voxels that belong to density rod itself and user-defined outermost voxels))
        for( size_t itr_radius( lower_radius_itr); itr_radius < upper_radius_itr; ++itr_radius)
        {
          // iterate over height in resolution steps
          for( size_t itr_length( 0); itr_length < m_Data.NumberLayers(); ++itr_length)
          {
            sum_intensities += m_Data( itr_length, itr_radius, itr_angle);
          }

          // initialize vector of correct height coordinates to pushback into 1D profile
          const double value_for_profile
          (
            ( itr_angle + 0.5) * m_AngleResolution
          );
          // pushback appropriate sum of intensities over the angles into 1D profile
          one_d_profile.PushBack( value_for_profile, sum_intensities);
        }
      }

      //return 1D profile
      return one_d_profile;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read MapCylindrical from std::istream
    //! @param ISTREAM std::istream from which the density is read
    //! @return std::istream
    std::istream &MapCylindrical::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_HeightResolution, ISTREAM); // height resolution value
      io::Serialize::Read( m_RadiusResolution, ISTREAM); // radius resolution value
      io::Serialize::Read( m_AngleResolution,  ISTREAM); // angle resolution value
      io::Serialize::Read( m_Minimum,          ISTREAM); // Minimal density value
      io::Serialize::Read( m_Maximum,          ISTREAM); // Maximal density value
      io::Serialize::Read( m_Mean,             ISTREAM); // Mean density value
      io::Serialize::Read( m_Body,             ISTREAM); // orientation of density as a body
      io::Serialize::Read( m_Data,             ISTREAM); // density as a tensor

      // end
      return ISTREAM;
    }

    //! @brief write MapCylindrical to std::ostream
    //! @param OSTREAM std::ostream to which the density is written
    //! @param INDENT indent used when writing the density map
    //! @return std::ostream
    std::ostream &MapCylindrical::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_HeightResolution,   OSTREAM, INDENT) << '\n'; // height resolution value
      io::Serialize::Write( m_RadiusResolution,   OSTREAM, INDENT) << '\n'; // radius resolution value
      io::Serialize::Write( m_AngleResolution,    OSTREAM, INDENT) << '\n'; // angle resolution value
      io::Serialize::Write( m_Minimum,            OSTREAM, INDENT) << '\n'; // Minimal density value
      io::Serialize::Write( m_Maximum,            OSTREAM, INDENT) << '\n'; // Maximal density value
      io::Serialize::Write( m_Mean,               OSTREAM, INDENT) << '\n'; // Mean density value

      // orientation of density as a body
      io::Serialize::Write( m_Body,      OSTREAM, INDENT) << '\n';

      // density as a tensor does also contain the dimensions of the density map
      io::Serialize::Write( m_Data,      OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace density
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
#include "density/bcl_density_mask_3d.h"

// includes from bcl - sorted alphabetically
#include "density/bcl_density_map.h"
#include "math/bcl_math_running_average_sd.h"
#include "util/bcl_util_si_ptr_vector.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Mask3d::s_Instance
    (
      GetObjectInstances().AddInstance( new Mask3d())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Mask3d::Mask3d() :
      m_Mask(),
      m_Index(),
      m_Position(),
      m_GridSpacing()
    {
    }

    //! @brief construct from list of coordinates
    //! @param COORDS a list of coordinates that defines the mask
    //! @param MASKING_DISTANCE distance for sigmoid in angstrom function to have marginal impact
    //! @param GRID_SPACING length of cube in grid
    //! @param CELL_POSITION position of the referenced grid in space
    Mask3d::Mask3d
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDS,
      const double MASKING_DISTANCE,
      const linal::Vector3D &GRID_SPACING,
      const linal::Vector3D &CELL_POSITION
    ) :
      m_Mask( 0, 0, 0),
      m_Index( 0, 0, 0),
      m_Position( CELL_POSITION),
      m_GridSpacing( GRID_SPACING)
    {
      // determine the corners as defined by given coordinates
      storage::VectorND< 2, linal::Vector3D> corners( DetermineGridCorners( COORDS));

      // undefined corners
      if( !corners.First().IsDefined())
      {
        BCL_MessageDbg( "could not determine corners for mask");
        return;
      }

      // add margin
      corners = AddMargin( corners, MASKING_DISTANCE);

      // find index of mask relative to map
      m_Index.First()  = int( std::floor( ( corners.First().X() - m_Position.X()) / GRID_SPACING.X()));
      m_Index.Second() = int( std::floor( ( corners.First().Y() - m_Position.Y()) / GRID_SPACING.Y()));
      m_Index.Third()  = int( std::floor( ( corners.First().Z() - m_Position.Z()) / GRID_SPACING.Z()));

      // allocate mask size and set values to 1
      m_Mask = math::Tensor< double>
      (
        size_t( std::ceil( ( corners.Second().X() - m_Position.X()) / GRID_SPACING.X())) - m_Index.First() ,
        size_t( std::ceil( ( corners.Second().Y() - m_Position.Y()) / GRID_SPACING.Y())) - m_Index.Second(),
        size_t( std::ceil( ( corners.Second().Z() - m_Position.Z()) / GRID_SPACING.Z())) - m_Index.Third() ,
        double( 1.0)
      );

      // store the real space index for easier access
      const linal::Vector3D realspaceindex
        (
          m_Position.X() + m_Index.First() * m_GridSpacing.X(),
          m_Position.Y() + m_Index.Second() * m_GridSpacing.Y(),
          m_Position.Z() + m_Index.Third() * m_GridSpacing.Z()
        );

      linal::Vector3D pos_voxel;

      // iterate over every point in grid to calculate the masking value
      for( size_t i( 0); i < m_Mask.GetNumberCols(); ++i)
      {
        pos_voxel.X() = i * m_GridSpacing.X() + realspaceindex.X();
        for( size_t j( 0); j < m_Mask.GetNumberRows(); ++j)
        {
          pos_voxel.Y() = j * m_GridSpacing.Y() + realspaceindex.Y();
          for( size_t k( 0); k < m_Mask.NumberLayers(); ++k)
          {
            pos_voxel.Z() = k * m_GridSpacing.Z() + realspaceindex.Z();

            double product( 1.0);
            // iterate over all coordinates
            for
            (
              util::SiPtrVector< const linal::Vector3D>::const_iterator itr( COORDS.Begin()), itr_end( COORDS.End());
              itr != itr_end;
              ++itr
            )
            {
              const linal::Vector3D &current_coord( **itr);
              if( !current_coord.IsDefined())
              {
                continue;
              }

              // calculate the distance between point in mask grid and specified coordinate
              product *= 1 - Sigmoid( MASKING_DISTANCE - ( current_coord - pos_voxel).Norm());
            }

            // actual mask weight is 1 minus the calculated product
            m_Mask( k, j, i) -= product;
          }
        }
      }
    };

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the cross correlation between two density maps weighted by mask
    //! @param DENSITY_MAP_EXP experimental density map
    //! @param DENSITY_MAP_SIM simulated density map
    //! @return CCC
    double Mask3d::CrossCorrelationCoefficient
    (
      const Map &DENSITY_MAP_EXP,
      const Map &DENSITY_MAP_SIM
    ) const
    {
      // assert the all maps' cell widths are identical
      BCL_Assert
      (
           DENSITY_MAP_EXP.GetCellWidth() == DENSITY_MAP_SIM.GetCellWidth()
        && DENSITY_MAP_EXP.GetCellWidth() == m_GridSpacing,
        "The density maps and the mask need to have the same Grid - spacing:/nmask: "
        + util::Format()(m_GridSpacing)
        + "\nsimulated density map: "
        + util::Format()( DENSITY_MAP_SIM.GetCellWidth())
        + "\nexperimental density map: "
        + util::Format()( DENSITY_MAP_EXP.GetCellWidth())
      );

      // determine overlap of each density map with this mask
      const storage::VectorND< 3, int> dm_rel_index_sim( RelativeIndex( DENSITY_MAP_SIM));
      const storage::VectorND< 3, int> dm_rel_index_exp( RelativeIndex( DENSITY_MAP_EXP));

      // determine the common overlap of the two densities and this mask
      const storage::VectorND< 2, storage::VectorND< 3, size_t> > commonarea
      (
        CommonOverlap( OverlappingIndices( DENSITY_MAP_EXP), OverlappingIndices( DENSITY_MAP_SIM))
      );

      // mean and sd of experimental and simulated map
      const math::RunningAverageSD< double> meansd_exp( CalculateMeanSDCommonRegion( DENSITY_MAP_EXP, commonarea));
      const math::RunningAverageSD< double> meansd_sim( CalculateMeanSDCommonRegion( DENSITY_MAP_SIM, commonarea));
      const double mean_exp( meansd_exp.GetAverage());
      const double mean_sim( meansd_sim.GetAverage());
      double coefficient( 0);

      // iterate over mask dimensions
      for( size_t i( commonarea.First().First()); i < commonarea.Second().First(); ++i)
      {
        const size_t exp_index_x( i - dm_rel_index_exp.First());
        const size_t sim_index_x( i - dm_rel_index_sim.First());
        for( size_t j( commonarea.First().Second()); j < commonarea.Second().Second(); ++j)
        {
          const size_t exp_index_y( j - dm_rel_index_exp.Second());
          const size_t sim_index_y( j - dm_rel_index_sim.Second());
          for( size_t k( commonarea.First().Third()); k < commonarea.Second().Third(); ++k)
          {
            const size_t exp_index_z( k - dm_rel_index_exp.Third());
            const size_t sim_index_z( k - dm_rel_index_sim.Third());
            // sum up products
            coefficient += m_Mask( k, j, i)
              *
              (
                DENSITY_MAP_EXP
                (
                  exp_index_x,
                  exp_index_y,
                  exp_index_z
                ) - mean_exp
              )
              *
              (
                DENSITY_MAP_SIM
                (
                  sim_index_x,
                  sim_index_y,
                  sim_index_z
                ) - mean_sim
              );
          }
        }
      }

      // number of voxels considered
      const size_t voxel_count
                (
                    ( commonarea.Second().First() - commonarea.First().First())
                  * ( commonarea.Second().Second() - commonarea.First().Second())
                  * ( commonarea.Second().Third() - commonarea.First().Third())
                );

      // calculate final correlation and return
      return coefficient * ( double( 1.0) / ( meansd_exp.GetStandardDeviation() * meansd_sim.GetStandardDeviation())) /
             voxel_count;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Mask3d::Read( std::istream &ISTREAM)
    {
      // read Index, Position and GridSpacing of the mask
      io::Serialize::Read( m_Index,        ISTREAM);
      io::Serialize::Read( m_Position,     ISTREAM);
      io::Serialize::Read( m_GridSpacing,  ISTREAM);

      // read the actual mask
      io::Serialize::Read( m_Mask,         ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &Mask3d::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write Index, Position and GridSpacing of the mask
      io::Serialize::Write( m_Index,        OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Position,     OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_GridSpacing,  OSTREAM, INDENT) << '\n';

      // write the actual mask
      io::Serialize::Write( m_Mask,         OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate mean and standard deviation of over the mask given a density map
    //! is not using the values in the mask, just the box that is defined by the mask
    //! @param DENSITY_MAP density map which's mean and SD shall be calculated for the area covered by the mask
    //! @return dataset statistics mean sd
    math::RunningAverageSD< double> Mask3d::CalculateMeanSD( const Map &DENSITY_MAP) const
    {
      // create empty data set statistic
      math::RunningAverageSD< double> mean_sd;

      // indices range
      const storage::VectorND< 2, storage::VectorND< 3, size_t> > beg_end_indices( OverlappingIndices( DENSITY_MAP));

      // relative index of density map
      const storage::VectorND< 3, int> dm_rel_index( RelativeIndex( DENSITY_MAP));

      // iterate over given map over mask boundaries
      for( size_t i( beg_end_indices.First().First()); i < ( beg_end_indices.Second().First()); ++i)
      {
        const size_t index_x( i - dm_rel_index.First());
        for( size_t j( beg_end_indices.First().Second()); j < beg_end_indices.Second().Second(); ++j)
        {
          const size_t index_y( j - dm_rel_index.Second());
          for( size_t k( beg_end_indices.First().Third()); k < beg_end_indices.Second().Third(); ++k)
          {
            const size_t index_z( k - dm_rel_index.Third());
            // fill mean_sd with values
            mean_sd += DENSITY_MAP( index_x, index_y, index_z);
          }
        }
      }

      // end
      return mean_sd;
    }

    //! @brief determine corners of grid the coordinates fit in
    //! @param COORDS a list of coordinates that defines the grid
    //! @return min coord xyz and max coord xyz
    storage::VectorND< 2, linal::Vector3D> Mask3d::DetermineGridCorners
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDS
    )
    {
      // Initialize corners
      const linal::Vector3D nullcoord( std::numeric_limits< double>::max());
      storage::VectorND< 2, linal::Vector3D> corners( nullcoord, -nullcoord);

      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator itr( COORDS.Begin()), itr_end( COORDS.End());
        itr != itr_end;
        ++itr
      )
      {
        const linal::Vector3D &current_coord( **itr);
        if( !current_coord.IsDefined())
        {
          // Make sure, that the current coordinates are defined
          BCL_MessageDbg
          (
            "cannot use undefined coordinates" + util::Format()( current_coord)
          );
          continue;
        }

        // iterate over all atoms to find the maximum and minimum values
        corners.First().X() = std::min( current_coord.X(), corners.First().X());
        corners.First().Y() = std::min( current_coord.Y(), corners.First().Y());
        corners.First().Z() = std::min( current_coord.Z(), corners.First().Z());

        corners.Second().X() = std::max( current_coord.X(), corners.Second().X());
        corners.Second().Y() = std::max( current_coord.Y(), corners.Second().Y());
        corners.Second().Z() = std::max( current_coord.Z(), corners.Second().Z());
      }

      // no coord was defined
      if( corners.First().X() > corners.Second().X())
      {
        corners.First() = linal::Vector3D( util::GetUndefined< double>());
        corners.Second() = linal::Vector3D( util::GetUndefined< double>());
      }

      return corners;
    }

    //! @brief add margins to corners
    //! @param CORNERS min coord xyz and max coord xyz
    //! @param MARGIN margin to add
    //! @return corners with additional margin
    storage::VectorND< 2, linal::Vector3D> Mask3d::AddMargin
    (
      const storage::VectorND< 2, linal::Vector3D> &CORNERS,
      const double MARGIN
    )
    {
      storage::VectorND< 2, linal::Vector3D> result( CORNERS); //Variable where the result shall be saved in
      result.First() -= MARGIN;
      result.Second() += MARGIN;

      return result;
    }

    //! @brief determine the indices ranges that overlap with given density map, relative to this mask
    //! @param DENSITY_MAP overlapping density map
    //! @return pair of start indices and end indices
    storage::VectorND< 2, storage::VectorND< 3, size_t> > Mask3d::OverlappingIndices( const Map &DENSITY_MAP) const
    {
      // relative index of density map
      const storage::VectorND< 3, int> dm_rel_index( RelativeIndex( DENSITY_MAP));

      // dimenesion of density map
      const storage::VectorND< 3, size_t> dm_dim( DENSITY_MAP.GetDimensions());

      // find start index iteration
      const int x_beg( std::max( dm_rel_index.First() , int( 0)));
      const int y_beg( std::max( dm_rel_index.Second(), int( 0)));
      const int z_beg( std::max( dm_rel_index.Third() , int( 0)));

      // find maximum values for each coordinate
      const int x_end( std::min( dm_rel_index.First()  + int( dm_dim.First()) , int( m_Mask.GetNumberCols())));
      const int y_end( std::min( dm_rel_index.Second() + int( dm_dim.Second()), int( m_Mask.GetNumberRows())));
      const int z_end( std::min( dm_rel_index.Third()  + int( dm_dim.Third()) , int( m_Mask.NumberLayers())));

      // if mask does not overlap with given density map
      if( x_beg >= x_end || y_beg >= y_end || z_beg >= z_end)
      {
        return storage::VectorND< 2, storage::VectorND< 3, size_t> >
               (
                 storage::VectorND< 3, size_t>( 0, 0, 0), storage::VectorND< 3, size_t>( 0, 0, 0)
               );
      }
      else
      {
        return storage::VectorND< 2, storage::VectorND< 3, size_t> >
               (
                 storage::VectorND< 3, size_t>( size_t( x_beg), size_t( y_beg), size_t( z_beg)),
                 storage::VectorND< 3, size_t>( size_t( x_end), size_t( y_end), size_t( z_end))
               );
      }
    }

    //! @brief calculate relative index of given density map to the mask
    //! @param DENSITY_MAP in question
    //! @return relative index
    storage::VectorND< 3, int> Mask3d::RelativeIndex( const Map &DENSITY_MAP) const
    {
      // index of density map
      // return
      return storage::VectorND< 3, int>
             (
               DENSITY_MAP.GetIndex()( 0) - m_Index.First(),
               DENSITY_MAP.GetIndex()( 1) - m_Index.Second(),
               DENSITY_MAP.GetIndex()( 2) - m_Index.Third()
             );
    }

    math::RunningAverageSD< double> Mask3d::CalculateMeanSDCommonRegion
    (
      const Map &DENSITY_MAP,
      const storage::VectorND< 2, storage::VectorND< 3, size_t> > &OVERLAP
    ) const
    {
      // create empty dataset
      math::RunningAverageSD< double> mean_sd;

      // relative index of density map
      const storage::VectorND< 3, int> dm_rel_index( RelativeIndex( DENSITY_MAP));

      for( size_t i( OVERLAP.First().First()); i < ( OVERLAP.Second().First()); ++i)
      {
        const size_t pos_x( i - dm_rel_index.First());
        for( size_t j( OVERLAP.First().Second()); j < OVERLAP.Second().Second(); ++j)
        {
          const size_t pos_y( j - dm_rel_index.Second());
          for( size_t k( OVERLAP.First().Third()); k < OVERLAP.Second().Third(); ++k)
          {
            const size_t pos_z( k - dm_rel_index.Third());
            // fill mean_sd with values
            mean_sd += DENSITY_MAP( pos_x, pos_y, pos_z);
          }
        }
      }
      return mean_sd;
    }

    //! @brief determines the overlap between the mask and two given ranges
    //! @param RANGE_ONE
    //! @param RANGE_TWO
    //! @return returns the common overlap
    storage::VectorND< 2, storage::VectorND< 3, size_t> > Mask3d::CommonOverlap
    (
      const storage::VectorND< 2, storage::VectorND< 3, size_t> > &RANGE_ONE,
      const storage::VectorND< 2, storage::VectorND< 3, size_t> > &RANGE_TWO
    )
    {
      storage::VectorND< 2, storage::VectorND< 3, size_t> > overlap;
      overlap.First().First()  = std::max( RANGE_ONE.First().First() , RANGE_TWO.First().First()) ;
      overlap.First().Second() = std::max( RANGE_ONE.First().Second(), RANGE_TWO.First().Second());
      overlap.First().Third()  = std::max( RANGE_ONE.First().Third() , RANGE_TWO.First().Third()) ;

      overlap.Second().First()  = std::min( RANGE_ONE.Second().First() , RANGE_TWO.Second().First()) ;
      overlap.Second().Second() = std::min( RANGE_ONE.Second().Second(), RANGE_TWO.Second().Second());
      overlap.Second().Third()  = std::min( RANGE_ONE.Second().Third() , RANGE_TWO.Second().Third()) ;

      return overlap;
    }

  } // namespace density
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
#include "density/bcl_density_protein_agreement_ccc.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinAgreementCCC::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinAgreementCCC( false, false))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param ADD_SIDECHAIN_ATOMS protein model will get side chains before the density simulation
    //! @param MULTIPLY_WITH_NUMBER_AAS multiply with number AAs so that the agreement scales with protein size
    ProteinAgreementCCC::ProteinAgreementCCC
    (
      const bool ADD_SIDECHAIN_ATOMS,
      const bool MULTIPLY_WITH_NUMBER_AAS
    ) :
      m_SimulatedContourLevelCutoff( 0.0),
      m_UseSideChains( ADD_SIDECHAIN_ATOMS),
      m_MultiplyWihtNumberAminoAcids( MULTIPLY_WITH_NUMBER_AAS)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new ProteinAgreementCCC copied from this one
    ProteinAgreementCCC *ProteinAgreementCCC::Clone() const
    {
      return new ProteinAgreementCCC( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinAgreementCCC::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinAgreementCCC::GetScheme() const
    {
      static const std::string s_scheme( "CCC");
      static const std::string s_scheme_sclaed( "CCCScaled");
      return m_MultiplyWihtNumberAminoAcids ? s_scheme_sclaed : s_scheme;
    }

    //! @brief access to the density simulator
    //! @return the density simulator used
    const util::ShPtr< SimulateInterface> &ProteinAgreementCCC::GetSimulator() const
    {
      return m_Simulate;
    }

    //! @brief set the simulator used for density agreement, e.g. for simulating density from the protein model
    //! @param SP_SIMULATOR ShPtr to SimulatInterface
    void ProteinAgreementCCC::SetSimulator( const util::ShPtr< SimulateInterface> &SP_SIMULATOR)
    {
      m_Simulate = SP_SIMULATOR;
    }

    //! @brief access to the density used for agreement calculation
    //! @return SiPtr to the density
    const util::SiPtr< const Map> &ProteinAgreementCCC::GetDensity() const
    {
      return m_Map;
    }

    //! @brief set the density used for agreement calculation
    //! @param SP_DENSITY SiPtr to the density map
    void ProteinAgreementCCC::SetDensityMap( const util::SiPtr< const Map> &SP_DENSITY)
    {
      m_Map = SP_DENSITY;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinAgreementCCC::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Calculates the standard deviation between a given protein model and a density map."
      );
      // serializer.AddInitializer
      // (
      //   "map",
      //   "density map to which the protein model is compared",
      //   io::Serialization::GetAgent( &m_Map)
      // );
      // serializer.AddInitializer
      // (
      //   "simulator",
      //   "algorithm to simulate a density map fora protein to facilitate comparison for CCC calculation.",
      //   io::Serialization::GetAgent( &m_Simulate)
      // );
      serializer.AddInitializer
      (
        "contour cutoff",
        "cutoff below which voxels are not considered for CCC calculation.",
        io::Serialization::GetAgent( &m_SimulatedContourLevelCutoff)
      );
      serializer.AddInitializer
      (
        "use side chains",
        "whether to add side chains to the protein model prior to CCC calculation.",
        io::Serialization::GetAgent( &m_UseSideChains)
      );
      serializer.AddInitializer
      (
        "multiply ccc",
        "multiply CCC with the number of amino acids in the protein model.",
        io::Serialization::GetAgent( &m_MultiplyWihtNumberAminoAcids)
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() taking the a ProteinModel and returning the cross correlation coefficient
    //! @param PROTEIN_MODEL
    //! @return correlation between the member density map and a simulated density map for PROTEIN_MODEL
    double ProteinAgreementCCC::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // density map
      Map simulated_map;

      if( m_UseSideChains)
      {
        // generate new protein model with side chains from PROTEINMODEL
        const util::ShPtr< assemble::ProteinModel> protein_model_with_side_chains
        (
          biol::AASideChainFactory( false, true).ProteinModelWithSideChains( PROTEIN_MODEL)
        );

        // simulate map
        simulated_map = m_Simulate->operator ()( protein_model_with_side_chains->GetAtoms());
      }
      else
      {
        // use protein as it is
        simulated_map = m_Simulate->operator ()( PROTEIN_MODEL.GetAtoms());
      }

      // cross correlation coefficient between simulated and actual map
      double ccc( CrossCorrelationCoefficient( *m_Map, simulated_map, m_SimulatedContourLevelCutoff));

      // scale to protein size
      if( m_MultiplyWihtNumberAminoAcids)
      {
        ccc *= PROTEIN_MODEL.GetNumberAAs();
      }

      // return correlation
      return -ccc;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &ProteinAgreementCCC::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Map                         , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Simulate                    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SimulatedContourLevelCutoff , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UseSideChains               , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MultiplyWihtNumberAminoAcids, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &ProteinAgreementCCC::Read( std::istream &ISTREAM)
    {
      // write member
      io::Serialize::Read( m_Map                         , ISTREAM);
      io::Serialize::Read( m_Simulate                    , ISTREAM);
      io::Serialize::Read( m_SimulatedContourLevelCutoff , ISTREAM);
      io::Serialize::Read( m_UseSideChains               , ISTREAM);
      io::Serialize::Read( m_MultiplyWihtNumberAminoAcids, ISTREAM);

      // end
      return ISTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the cross correlation between experimental and simulated density map
    //! this ccc measure only considers voxels, if the intensity in the simulated map is above the given contour level
    //! this prevents that adjacent protein densities in experimental maps have a negative contribution to the measure
    //! @param EXPERIMENTAL_DENSITY_MAP map from experiment
    //! @param SIMULATED_DENSITY_MAP map simulated from protein structure
    //! @param CONTOUR_LEVEL_SIMULATED the minimal intensity in the simulated map voxel, to be considered for the CCC
    //! @return cross correlation coefficient
    double ProteinAgreementCCC::CrossCorrelationCoefficient
    (
      const Map &EXPERIMENTAL_DENSITY_MAP,
      const Map &SIMULATED_DENSITY_MAP,
      const double CONTOUR_LEVEL_SIMULATED
    )
    {
      // create common sub tensor
      const storage::VectorND< 2, math::Tensor< double> > exp_sim_sub_tensor( EXPERIMENTAL_DENSITY_MAP.CommonSubTensor( SIMULATED_DENSITY_MAP));

      // number of voxels that are above the CONTOUR_LEVEL in the experimental (this) density map
      size_t count_voxel( 0);
      double sum_sim( 0);
      double sum_exp( 0);
      double sum_sim2( 0);
      double sum_exp2( 0);
      double sum_exp_sim( 0);

      // iterate over experimental and simulated sub tensor
      for
      (
        const double
          *exp( exp_sim_sub_tensor.First().Begin()), *exp_end( exp_sim_sub_tensor.First().End()),
          *sim( exp_sim_sub_tensor.Second().Begin());
        exp != exp_end;
        ++exp, ++sim
      )
      {
        const double sim_int( *sim);

        // ignore sim and exp intensities below sim contour level
        if( sim_int > CONTOUR_LEVEL_SIMULATED)
        {
          const double exp_int( *exp);
          ++count_voxel;
          sum_exp += exp_int;
          sum_sim += sim_int;
          sum_exp_sim += exp_int * sim_int;
          sum_exp2 += math::Sqr( exp_int);
          sum_sim2 += math::Sqr( sim_int);
        }
      }

      // calculate actual correlation
      double correlation( count_voxel * sum_exp_sim - sum_exp * sum_sim);
      correlation /= math::Sqrt( count_voxel * sum_exp2 - math::Sqr( sum_exp)) * math::Sqrt( count_voxel * sum_sim2 - math::Sqr( sum_sim));

      // end
      return correlation;
    }

  } // namespace density
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
#include "density/bcl_density_protein_agreement_likelihood.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_mask_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
  //////////
  // data //
  //////////

    //! @brief low res masking distance for Mask3D using CA only
    const double ProteinAgreementLikelihood::s_LowResolutionMaskingDistance( 8.0);

    //! @brief high res masking distance for Mask3D using all side chain atoms
    const double ProteinAgreementLikelihood::s_HighResolutionMaskingDistance( 5.0);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinAgreementLikelihood::s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinAgreementLikelihood::ProteinAgreementLikelihood() :
      m_LogLikelihood(),
      m_HighResolution( false),
      m_AtomTypes( storage::Set< biol::AtomType>( biol::GetAtomTypes().CA)),
      m_DensityMap(),
      m_Simulator()
    {
      SetScheme();
    }

    //! @brief constructor from atom types
    //! @param HIGH_RESOLUTION use all atoms of residue plus neighboring residues, otherwise just CA
    //! @param ATOM_TYPES atom types to consider
    //! @param MEAN_CCC_RESOLUTION_FUNCTION
    //! @param SD_CCC_RESOLUTION_FUNCTION
    ProteinAgreementLikelihood::ProteinAgreementLikelihood
    (
      const bool HIGH_RESOLUTION,
      const storage::Set< biol::AtomType> &ATOM_TYPES,
      const math::FunctionInterfaceSerializable< double, double> &MEAN_CCC_RESOLUTION_FUNCTION,
      const math::FunctionInterfaceSerializable< double, double> &SD_CCC_RESOLUTION_FUNCTION
    ) :
      m_HighResolution( HIGH_RESOLUTION),
      m_AtomTypes( ATOM_TYPES),
      m_MeanFit( MEAN_CCC_RESOLUTION_FUNCTION.Clone()),
      m_SdFit( SD_CCC_RESOLUTION_FUNCTION.Clone())
    {
      SetScheme();
    }

    //! @brief Clone function
    //! @return pointer to new ProteinAgreementLikelihood
    ProteinAgreementLikelihood *ProteinAgreementLikelihood::Clone() const
    {
      return new ProteinAgreementLikelihood( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinAgreementLikelihood::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinAgreementLikelihood::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief access to the density simulator
    //! @return the density simulator used
    const util::ShPtr< SimulateInterface> &ProteinAgreementLikelihood::GetSimulator() const
    {
      return m_Simulator;
    }

    //! @brief set the simulator used for density agreement, e.g. for simulating density from the protein model
    //! @param SP_SIMULATOR ShPtr to SimulatInterface
    void ProteinAgreementLikelihood::SetSimulator( const util::ShPtr< SimulateInterface> &SP_SIMULATOR)
    {
      m_Simulator = SP_SIMULATOR;
      const double resolution( SP_SIMULATOR->GetResolution());
      m_LogLikelihood = math::LogLikelihood( m_MeanFit->operator ()( resolution), m_SdFit->operator ()( resolution));
    }

    //! @brief access to the density used for agreement calculation
    //! @return SiPtr to the density
    const util::SiPtr< const Map> &ProteinAgreementLikelihood::GetDensity() const
    {
      return m_DensityMap;
    }

    //! @brief set the density used for agreement calculation
    //! @param SP_DENSITY SiPtr to the density map
    void ProteinAgreementLikelihood::SetDensityMap( const util::SiPtr< const Map> &SP_DENSITY)
    {
      m_DensityMap = SP_DENSITY;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the standard Deviation between the protein model and the given density
    //! @param PROTEIN_MODEL protein of interest
    //! @return score for the density likelihood
    double ProteinAgreementLikelihood::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // high resolution score
      if( m_HighResolution)
      {
        return HighResolutionScore( PROTEIN_MODEL);
      }
      // low resolution score
      else
      {
        return LowResolutionScore( PROTEIN_MODEL);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinAgreementLikelihood::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DensityMap    , ISTREAM);
      io::Serialize::Read( m_HighResolution, ISTREAM);
      io::Serialize::Read( m_AtomTypes     , ISTREAM);
      io::Serialize::Read( m_LogLikelihood , ISTREAM);
      io::Serialize::Read( m_Simulator     , ISTREAM);

      // set scheme according to atom types
      SetScheme();

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinAgreementLikelihood::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_DensityMap    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HighResolution, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomTypes     , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_LogLikelihood , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Simulator     , OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate low resolution score, only considering CA atoms
    //! @param PROTEIN_MODEL protein of interest
    //! @return score for the density likelihood
    double ProteinAgreementLikelihood::LowResolutionScore( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // simulate density for given PROTEIN_MODEL
      const Map density_sim( m_Simulator->operator()( PROTEIN_MODEL.GetAtoms( m_AtomTypes)));

      // get all residues
      const util::SiPtrVector< const biol::AABase> residues( PROTEIN_MODEL.GetAminoAcids());

      // score
      double score( 0.0);

      // iterate over all residues
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator aa_itr( residues.Begin()), aa_itr_end( residues.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // current res coordinates
        const util::SiPtrVector< const linal::Vector3D> coords( ( *aa_itr)->GetAtomCoordinates( m_AtomTypes));

        // skip undefined, as for GLY
        if( !coord::AreDefinedCoordinates( coords))
        {
          continue;
        }

        // generate mask with coords
        const Mask3d current_mask
        (
          coords, s_LowResolutionMaskingDistance, density_sim.GetCellWidth(), density_sim.GetOrigin()
        );

        // calculate cross correlation over mask
        const double mask_ccc( current_mask.CrossCorrelationCoefficient( *m_DensityMap, density_sim));
        if( util::IsDefined( mask_ccc))
        {
          // evaluate log likelihood function and add to score
          score += m_LogLikelihood( mask_ccc);
        }
      }

      // end
      return score;
    }

    //! @brief calculate high resolution score, considering all atoms of each amino acid sidechain plus the two
    //! neighboring residues
    //! @param PROTEIN_MODEL protein of interest
    //! @return score for the density likelihood
    double ProteinAgreementLikelihood::HighResolutionScore( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // generate new protein model with side chains from PROTEIN_MODEL
      const assemble::ProteinModel model
      (
        *biol::AASideChainFactory( false, true).ProteinModelWithSideChains( PROTEIN_MODEL)
      );

      // get all residues
      const util::SiPtrVector< const biol::AABase> residues( model.GetAminoAcids());

      // score
      double score( 0.0);

      // need at least 3 amino acids for high res score
      if( residues.GetSize() < 3)
      {
        return score;
      }

      // simulate density for given PROTEIN_MODEL
      const Map density_sim( m_Simulator->operator()( model.GetAtoms()));

      // iterate over all residues
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          aa_itr_left( residues.Begin()), aa_itr_current( aa_itr_left + 1), aa_itr_right( aa_itr_current + 1),
          aa_itr_end( residues.End());
        aa_itr_right != aa_itr_end;
        ++aa_itr_left, ++aa_itr_current, ++aa_itr_right
      )
      {
        const biol::AABase &aa_left( **aa_itr_left);
        const biol::AABase &aa_current( **aa_itr_current);
        const biol::AABase &aa_right( **aa_itr_right);

        // sequence separation
        const size_t seq_separation_lc( biol::SequenceSeparation( aa_left, aa_current));
        const size_t seq_separation_cr( biol::SequenceSeparation( aa_current, aa_right));

        // check that all amino acids are from the same chain and are neighboring amino acids
        if( seq_separation_lc != 0 || seq_separation_cr != 0)
        {
          continue;
        }

        util::SiPtrVector< const linal::Vector3D> coords;
        coords.Append( aa_left.GetAtomCoordinates());
        coords.Append( aa_current.GetAtomCoordinates());
        coords.Append( aa_right.GetAtomCoordinates());
        const Mask3d current_mask
        (
          coords, s_HighResolutionMaskingDistance, density_sim.GetCellWidth(), density_sim.GetOrigin()
        );

        // calculate cross correlation over mask
        const double mask_ccc( current_mask.CrossCorrelationCoefficient( *m_DensityMap, density_sim));
        if( util::IsDefined( mask_ccc))
        {
          // evaluate log likelihood function and add to score
          score += m_LogLikelihood( mask_ccc);
        }
      }

      // end
      return score;
    }

    //! @brief string for the default scheme
    const std::string &ProteinAgreementLikelihood::GetDefaultScheme() const
    {
      // initialize static scheme
      static const std::string s_scheme( "Likelihood");

      // end
      return s_scheme;
    }

    //! @brief set scheme
    //! @details set scheme from the default scheme and atom types
    void ProteinAgreementLikelihood::SetScheme()
    {
      m_Scheme = GetDefaultScheme();
      for
      (
        storage::Set< biol::AtomType>::const_iterator itr( m_AtomTypes.Begin()), itr_end( m_AtomTypes.End());
        itr != itr_end;
        ++itr
      )
      {
        m_Scheme += ( *itr)->GetName();
      }
    }

  } // namespace density
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
#include "density/bcl_density_protein_agreements.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom_types.h"
#include "density/bcl_density_protein_agreement_ccc.h"
#include "density/bcl_density_protein_agreement_likelihood.h"
#include "density/bcl_density_simulators.h"
#include "math/bcl_math_polynomial.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

    const double g_LikelihoodCaMeanPolynomCoeff[]       = {  0.00000346634, -0.000149827,   0.00159432, 0.00788547, 0.00386209};
    const double g_LikelihoodCaSdPolynomCoeff[]         = { -0.000000465365, 0.0000547865, -0.00153436, 0.0161681,  0.0339695};
    const double g_LikelihoodCbMeanPolynomCoeff[]       = {  0.00000322947, -0.000139397,   0.00145854, 0.00815733, 0.00711573};
    const double g_LikelihoodCbSdPolynomCoeff[]         = { -0.00000126792,  0.0000982248, -0.00234045, 0.0217593,  0.0257755};
    const double g_LikelihoodBackBoneMeanPolynomCoeff[] = {  0.00000427199,  0.000235853,  -0.00546931, 0.0659212, -0.12085};
    const double g_LikelihoodBackBoneSdPolynomCoeff[]   = { -0.00005495,     0.000288613,  -0.00544743, 0.044201,  -0.021971};

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinAgreements::ProteinAgreements() :
      e_CCC(                AddAgreement( util::ShPtr< ProteinAgreementInterface>( new ProteinAgreementCCC( false, false)))),
      e_CCCScaled(          AddAgreement( util::ShPtr< ProteinAgreementInterface>( new ProteinAgreementCCC( false, true)))),
      e_LikelihoodCa(       AddAgreement( util::ShPtr< ProteinAgreementInterface>( new ProteinAgreementLikelihood( false, storage::Set< biol::AtomType>( biol::GetAtomTypes().CA), math::Polynomial::MakeFromCoefficients( linal::Vector< double>( 5, g_LikelihoodCaMeanPolynomCoeff)), math::Polynomial::MakeFromCoefficients( linal::Vector< double>( 5, g_LikelihoodCaSdPolynomCoeff)))))),
      e_LikelihoodCb(       AddAgreement( util::ShPtr< ProteinAgreementInterface>( new ProteinAgreementLikelihood( false, storage::Set< biol::AtomType>( biol::GetAtomTypes().CB), math::Polynomial::MakeFromCoefficients( linal::Vector< double>( 5, g_LikelihoodCbMeanPolynomCoeff)), math::Polynomial::MakeFromCoefficients( linal::Vector< double>( 5, g_LikelihoodCbSdPolynomCoeff)))))),
      e_LikelihoodBackBone( AddAgreement( util::ShPtr< ProteinAgreementInterface>( new ProteinAgreementLikelihood( false, biol::GetAtomTypes().GetBackBoneAtomTypes()            , math::Polynomial::MakeFromCoefficients( linal::Vector< double>( 5, g_LikelihoodBackBoneMeanPolynomCoeff)), math::Polynomial::MakeFromCoefficients( linal::Vector< double>( 5, g_LikelihoodBackBoneSdPolynomCoeff))))))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinAgreements::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add a new protein agreement instance
    //! @param SP_PROTEIN_AGREEMENT shptr to instance of a ProteinAgreementInterface derived class
    //! @return the new enum constructed from the agreement
    ProteinAgreement ProteinAgreements::AddAgreement( const util::ShPtr< ProteinAgreementInterface> &SP_PROTEIN_AGREEMENT)
    {
      return AddEnum( SP_PROTEIN_AGREEMENT->GetScheme(), SP_PROTEIN_AGREEMENT);
    }

    //! @brief create ProteinAgreement from enum, density map and resolution
    //! @param PROTEIN_AGREEMENT enum of the agreement to be used
    //! @param SIMULATOR
    //! @param SP_DENSITY density map to be used
    //! @param RESOLUTION resolution in Angstrom to be used
    //! @return ShPtr to density::ProteinAgreement
    util::ShPtr< ProteinAgreementInterface> ProteinAgreements::CreateProteinAgreement
    (
      const ProteinAgreement &PROTEIN_AGREEMENT,
      const Simulator &SIMULATOR,
      const util::SiPtr< const Map> &SP_DENSITY,
      const double RESOLUTION
    ) const
    {
      // check that the enum is defined
      if( !PROTEIN_AGREEMENT.IsDefined() || !PROTEIN_AGREEMENT->IsDefined())
      {
        return util::ShPtr< ProteinAgreementInterface>();
      }

      // create simulator
      util::ShPtr< SimulateInterface> sp_simulator( GetSimulators().CreateSimulator( SIMULATOR, SP_DENSITY->GetCellWidth(), RESOLUTION));

      // clone an agreement
      util::ShPtr< ProteinAgreementInterface> sp_protein_agreement( PROTEIN_AGREEMENT->HardCopy());

      // set the member
      sp_protein_agreement->SetDensityMap( SP_DENSITY);
      sp_protein_agreement->SetSimulator( sp_simulator);

      // end
      return sp_protein_agreement;
    }

    //! @brief construct on access function for all ProteinAgreements
    //! @return reference to only instances of ProteinAgreements
    ProteinAgreements &GetProteinAgreements()
    {
      return ProteinAgreements::GetEnums();
    }

  } // namespace density

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< density::ProteinAgreementInterface>, density::ProteinAgreements>;

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
#include "density/bcl_density_simulate_default.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SimulateDefault::s_Instance
    (
      GetObjectInstances().AddInstance( new SimulateDefault( linal::Vector3D( 2.3, 2.3, 2.3), 6.9, e_Gaussian))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from grid spacing and resolution
    //! @param GRID_SPACING the spacing for the density grid
    //! @param RESOLUTION the resolution to simulate for
    //! @param SMOOTHING_KERNEL kernel for the smoothing
    SimulateDefault::SimulateDefault
    (
      const linal::Vector3D &GRID_SPACING,
      const double RESOLUTION,
      const Kernel SMOOTHING_KERNEL
    ) :
      m_GridSpacing( GRID_SPACING),
      m_Resolution( RESOLUTION),
      m_Margin( 2),
      m_SmoothingKernel( SMOOTHING_KERNEL)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SimulateDefault
    SimulateDefault *SimulateDefault::Clone() const
    {
      return new SimulateDefault( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SimulateDefault::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the resolution
    //! @param RESOLUTION the resolution for the density map to be generated
    void SimulateDefault::SetResolution( const double RESOLUTION)
    {
      m_Resolution = RESOLUTION;
    }

    //! @brief set the resolution
    double SimulateDefault::GetResolution() const
    {
      return m_Resolution;
    }

    //! @brief set the grid spacing
    //! @param GRID_SPACING the width of a grid element in x, y and z
    void SimulateDefault::SetGridSpacing( const linal::Vector3D &GRID_SPACING)
    {
      m_GridSpacing = GRID_SPACING;
    }

    //! @brief set the margin
    //! @param MARGIN number of additional cells next to last atom occupied cells
    void SimulateDefault::SetMargin( const size_t MARGIN)
    {
      m_Margin = MARGIN;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief generate simulated density from given list of atoms
    //! @param ATOMS siptrvector of atoms
    //! @return a simulated density map
    Map SimulateDefault::operator()( const util::SiPtrVector< const biol::Atom> &ATOMS) const
    {
      // check that any atoms are given
      if( ATOMS.IsEmpty())
      {
        return Map();
      }

      // measure ATOMS extent
      storage::VectorND< 2, linal::Vector3D> min_max_coord( DetermineGridCorners( ATOMS));

      linal::Vector3D &mincoord( min_max_coord.First());
      linal::Vector3D &maxcoord( min_max_coord.Second());

      // bring lattice into register with origin
      for( size_t i( 0); i < 3; ++i)
      {
        mincoord( i) = m_GridSpacing( i) * std::floor( mincoord( i) / m_GridSpacing( i));
        maxcoord( i) = m_GridSpacing( i) * std::ceil( maxcoord( i) / m_GridSpacing( i));
      }

      // allocate atom density map
      const storage::VectorND< 3, size_t> ext_atom_map // extension of atom map
      (
        size_t( ceil( ( maxcoord.X() - mincoord.X()) / m_GridSpacing.X()) + 2 * m_Margin + 1),
        size_t( ceil( ( maxcoord.Y() - mincoord.Y()) / m_GridSpacing.Y()) + 2 * m_Margin + 1),
        size_t( ceil( ( maxcoord.Z() - mincoord.Z()) / m_GridSpacing.Z()) + 2 * m_Margin + 1)
      );
      const linal::Vector3D grid2
      (
        mincoord.X() - m_Margin * m_GridSpacing.X(),
        mincoord.Y() - m_Margin * m_GridSpacing.Y(),
        mincoord.Z() - m_Margin * m_GridSpacing.Z()
      );

      math::Tensor< double> atom_map( ext_atom_map.Third(), ext_atom_map.Second(), ext_atom_map.First(), double( 0));

      // correct for lattice interpolation smoothing effects
      // slightly lowers the kernel width to maintain target resolution
      const size_t corrmode( 2);

      // desired kernel amplitude (scaling factor):
      const double kernel_amplitude( 1.0);

      // interpolate structure to protein map and keep track of variability
      // Projecting atoms to cubic lattice by trilinear interpolation...
      double varp( 0.0);
      size_t total_weight( 0);
      for
      (
        util::SiPtrVector< const biol::Atom>::const_iterator
          atom_itr( ATOMS.Begin()), atom_itr_end( ATOMS.End());
        atom_itr != atom_itr_end;
        ++atom_itr
      )
      {
        const biol::Atom &current_atom( **atom_itr);
        const linal::Vector3D &current_coord( current_atom.GetCoordinates());

        if( !current_coord.IsDefined())
        {
          continue;
        }
        // weigh by atom type
        size_t weight( 0);
        if( !current_atom.GetType()->GetElementType().IsDefined())
        {
          continue;
        }
        else if( current_atom.GetType()->GetElementType() <= chemistry::GetElementTypes().e_Lithium)  continue;   // till Lithium
        else if( current_atom.GetType()->GetElementType() <= chemistry::GetElementTypes().e_Neon)     weight = 1; // till Neon
        else if( current_atom.GetType()->GetElementType() <= chemistry::GetElementTypes().e_Sulfur)   weight = 2; // till Sulfur
        else if( current_atom.GetType()->GetElementType() <= chemistry::GetElementTypes().e_Titanium) weight = 3; // till Titanium
        else if( current_atom.GetType()->GetElementType() <= chemistry::GetElementTypes().e_Nickel)   weight = 4; // till Nickel
        else weight = 5;                                             // from Copper

        // sum up weights
        total_weight += weight;

        // compute position within grid in Angstroem
        const linal::Vector3D grid_position
                             (
                               m_Margin + ( current_coord.X() - mincoord.X()) / m_GridSpacing.X(),
                               m_Margin + ( current_coord.Y() - mincoord.Y()) / m_GridSpacing.Y(),
                               m_Margin + ( current_coord.Z() - mincoord.Z()) / m_GridSpacing.Z()
                             );
        // convert to
        const int x0( int( std::floor( grid_position.X())));
        const int y0( int( std::floor( grid_position.Y())));
        const int z0( int( std::floor( grid_position.Z())));

        // interpolate - determine position in the box that point is interpolated to
        const double a( x0 + 1 - grid_position.X());
        const double b( y0 + 1 - grid_position.Y());
        const double c( z0 + 1 - grid_position.Z());

        const double ab( a * b);
        const double ac( a * c);
        const double bc( b * c);
        const double aa( a * a);
        const double bb( b * b);
        const double cc( c * c);
        const double oma( 1 - a);
        const double omb( 1 - b);
        const double omc( 1 - c);
        const double omab( oma * omb);
        const double omac( oma * omc);
        const double ombc( omb * omc);
        const double omaa( oma * oma);
        const double ombb( omb * omb);
        const double omcc( omc * omc);
        const double val1(   ab *    c * weight);
        const double val2(   ab *  omc * weight);
        const double val3(   ac *  omb * weight);
        const double val4(   bc *  oma * weight);
        const double val5(    a * ombc * weight);
        const double val6(    c * omab * weight);
        const double val7(    b * omac * weight);
        const double val8( omab * omc  * weight);

        atom_map( z0    , y0    , x0    ) += val1;
        atom_map( z0 + 1, y0    , x0    ) += val2;
        atom_map( z0    , y0 + 1, x0    ) += val3;
        atom_map( z0    , y0    , x0 + 1) += val4;
        atom_map( z0 + 1, y0 + 1, x0    ) += val5;
        atom_map( z0    , y0 + 1, x0 + 1) += val6;
        atom_map( z0 + 1, y0    , x0 + 1) += val7;
        atom_map( z0 + 1, y0 + 1, x0 + 1) += val8;

        varp += val1 * ( omaa + ombb + omcc);
        varp += val2 * ( omaa + ombb +   cc);
        varp += val3 * ( omaa +   bb + omcc);
        varp += val4 * (   aa + ombb + omcc);
        varp += val5 * ( omaa +   bb +   cc);
        varp += val6 * (   aa +   bb + omcc);
        varp += val7 * (   aa + ombb +   cc);
        varp += val8 * (   aa +   bb +   cc);
      }

      // normalize
      varp /= double( total_weight);

      // target resolution (2 sigma)
      const double sigma( m_Resolution / 2);

      linal::Vector< double> unknown_fac( s_MaxKernel, double( 0));
      unknown_fac( e_Gaussian)          =  0.0;
      unknown_fac( e_Triangular)        =  1.0;
      unknown_fac( e_SemiEpanechnikov)  =  1.5;
      unknown_fac( e_Epanechnikov)      =  2.0;
      unknown_fac( e_HardSphere)        = 60.0;

      linal::Vector< double> radius_half( s_MaxKernel, double( 0));
      radius_half( e_Gaussian)         = sigma * sqrt( log( 2.0)) / sqrt(1.5);
      radius_half( e_Triangular)       = sigma / ( exp( ( 1.0 / unknown_fac( e_Triangular))       * log( 2.0)) * sqrt( 3.0*( 3.0+unknown_fac( e_Triangular))       / (5.0*(5.0+unknown_fac( e_Triangular)))));
      radius_half( e_SemiEpanechnikov) = sigma / ( exp( ( 1.0 / unknown_fac( e_SemiEpanechnikov)) * log( 2.0)) * sqrt( 3.0*( 3.0+unknown_fac( e_SemiEpanechnikov)) / (5.0*(5.0+unknown_fac( e_SemiEpanechnikov)))));
      radius_half( e_Epanechnikov)     = sigma / ( exp( ( 1.0 / unknown_fac( e_Epanechnikov))     * log( 2.0)) * sqrt( 3.0*( 3.0+unknown_fac( e_Epanechnikov))     / (5.0*(5.0+unknown_fac( e_Epanechnikov)))));
      radius_half( e_HardSphere)       = sigma / ( exp( ( 1.0 / unknown_fac( e_HardSphere))       * log( 2.0)) * sqrt( 3.0*( 3.0+unknown_fac( e_HardSphere))       / (5.0*(5.0+unknown_fac( e_HardSphere)))));

      linal::Vector< double> radius_cut( s_MaxKernel, double( 0));
      radius_cut( e_Gaussian)         = sqrt( 3.0) * sigma;
      radius_cut( e_Triangular)       = ( exp( ( 1.0 / unknown_fac( e_Triangular))       * log( 2.0))) * radius_half( e_Triangular);
      radius_cut( e_SemiEpanechnikov) = ( exp( ( 1.0 / unknown_fac( e_SemiEpanechnikov)) * log( 2.0))) * radius_half( e_SemiEpanechnikov);
      radius_cut( e_Epanechnikov)     = ( exp( ( 1.0 / unknown_fac( e_Epanechnikov))     * log( 2.0))) * radius_half( e_Epanechnikov);
      radius_cut( e_HardSphere)       = ( exp( ( 1.0 / unknown_fac( e_HardSphere))       * log( 2.0))) * radius_half( e_HardSphere);

      Kernel kernel( m_SmoothingKernel);

      // if the half radius is too small than apply no smoothing
      if( radius_half( e_Gaussian) / m_GridSpacing.X() < 1.0)
      {
        kernel = e_NoSmoothing;
      }

      double kmsd( 0);
      if( kernel == e_Gaussian)
      {
        kmsd = math::Sqr( sigma) / ( m_GridSpacing.X() * m_GridSpacing.X());
      }
      else
      {
        kmsd = math::Sqr( sigma) / math::Sqr( m_GridSpacing.X());
      }

      double varmap( 0);
      if( corrmode == 1)
      {
        varmap -= varp; // variances are additive for uncorrelated samples
      }
      else // ( corrmode == 2)
      {
        varmap = kmsd;
      }

      if( varmap < 0)
      {
        // lattice smoothing exceeds kernel size
        return Map();
      }

      // maximal spanning of density coming from one atom
      int exth( 0);
      double sigmamap( 0);

      // compute lattice noise corrected kernel maps
      switch( kernel)
      {
        case e_NoSmoothing:
          break;

        case e_Gaussian:
        {
          sigmamap = sqrt( varmap / 3.0);   // sigma-1D
          exth = int( ceil( 3 * sigmamap)); // truncate at 3 sigma-1D == sqrt(3) sigma-3D
          break;
        }

        case e_Triangular:
        case e_SemiEpanechnikov:
        case e_Epanechnikov:
        case e_HardSphere:
        {
          exth = int( ceil( radius_cut( kernel) / m_GridSpacing.X()));
          break;
        }

        default:
          break;
      }

      // allocate kernel map
      math::Tensor< double> kernel_map( 2 * exth + 1, 2 * exth + 1, 2 * exth + 1, double( 0));

      // smooth map depending on kernel
      switch( kernel)
      {
        case e_NoSmoothing:
          break;

        case e_Gaussian:
        {
          // write Gaussian within 3 sigma-1D to map
          const double bvalue( -1.0 / ( 2.0 * math::Sqr( sigmamap)));
          const double cvalue( 9 * math::Sqr( sigmamap));
          for( size_t indz( 0); indz < kernel_map.NumberLayers(); ++indz)
          {
            for( size_t indy( 0); indy < kernel_map.GetNumberRows(); ++indy)
            {
              for( size_t indx( 0); indx < kernel_map.GetNumberCols(); ++indx)
              {
                const double dsqu( math::Sqr( indx - exth) + math::Sqr( indy - exth) + math::Sqr( indz - exth));
                if( dsqu < cvalue)
                {
                  kernel_map( indz, indy, indx) = kernel_amplitude * exp( dsqu * bvalue);
                }
              }
            }
          }

          break;
        }

        case e_Triangular:
        case e_SemiEpanechnikov:
        case e_Epanechnikov:
        case e_HardSphere:
        {
          // write kernel to map
          const double bvalue( 0.5 * exp( -unknown_fac( kernel) * log( radius_half( kernel))) * exp( unknown_fac( kernel) * log( m_GridSpacing.X())));
          for( size_t indz( 0); indz < kernel_map.NumberLayers(); ++indz)
          {
            for( size_t indy( 0); indy < kernel_map.GetNumberRows(); ++indy)
            {
              for( size_t indx( 0); indx < kernel_map.GetNumberCols(); ++indx)
              {
                const double dsqu
                (
                  exp
                  (
                    ( unknown_fac( kernel) / 2.0)
                    * std::log
                    (
                      double
                      (
                          math::Sqr( indx - exth)
                        + math::Sqr( indy - exth)
                        + math::Sqr( indz - exth)
                      )
                    )
                  )
                );
                kernel_map( indz, indy, indx) = std::max( double( 0.0), kernel_amplitude * ( 1.0 - dsqu * bvalue));
              }
            }
          }
          break;
        }

        default:
          break;
      }

      // convolve and write output
      switch( kernel)
      {
        case e_NoSmoothing:
        {
          // index of map
          const linal::VectorND< int, 3> index
          (
            int( floor( grid2.X() / m_GridSpacing.X() + 0.5)),
            int( floor( grid2.Y() / m_GridSpacing.Y() + 0.5)),
            int( floor( grid2.Z() / m_GridSpacing.Z() + 0.5))
          );

          // intervals
          linal::VectorND< int, 3> intervals;
          for( size_t i( 0); i < 3; ++i)
          {
            if( index( i) <= -int( ext_atom_map( i)))
            {
              intervals( i) = math::Absolute( index( i));
            }
            else if( index( i) >= 0)
            {
              intervals( i) = index( i) + ext_atom_map( i) - 1;
            }
            else
            {
              intervals( i) = ext_atom_map( i) - 1;
            }
          }

          // length
          const linal::Vector3D length
          (
            m_GridSpacing.X() * double( intervals( 0)),
            m_GridSpacing.Y() * double( intervals( 1)),
            m_GridSpacing.Z() * double( intervals( 2))
          );

//          BCL_Message
//          (
//            util::Message::e_Standard,
//            util::Format()( "kernel width should be large relative to lattice voxel spacing!\n") +
//            "Spatial resolution (2 sigma) from lattice smoothing alone: " +
//            util::Format().W( 6).FFP( 3)( m_GridSpacing.X() * sqrt(varp) * 2.0) + "\n" +
//            "If kernel smoothing desired, increase kernel width or decrease lattice spacing."
//          );

          return Map
          (
            atom_map,
            index,
            intervals,
            length,
            m_GridSpacing,
            Map::GetDefaultAngle(),
            Map::GetDefaultAxis(),
            linal::Vector3D( 0.0)
          );
        }
        case e_Gaussian:
        case e_Triangular:
        case e_SemiEpanechnikov:
        case e_Epanechnikov:
        case e_HardSphere:
        {
          // Convolving lattice with kernel

          // allocate output density map
          const storage::VectorND< 3, size_t> ext_output // extension output map
          (
            size_t( ceil( ( maxcoord.X() - mincoord.X()) / m_GridSpacing.X()) + 2 * exth + 2 * m_Margin + 1),
            size_t( ceil( ( maxcoord.Y() - mincoord.Y()) / m_GridSpacing.Y()) + 2 * exth + 2 * m_Margin + 1),
            size_t( ceil( ( maxcoord.Z() - mincoord.Z()) / m_GridSpacing.Z()) + 2 * exth + 2 * m_Margin + 1)
          );
          const linal::Vector3D grid3
          (
            mincoord.X() - ( exth + m_Margin) * m_GridSpacing.X(),
            mincoord.Y() - ( exth + m_Margin) * m_GridSpacing.Y(),
            mincoord.Z() - ( exth + m_Margin) * m_GridSpacing.Z()
          );
          math::Tensor< double> output_map( ext_output.Third(), ext_output.Second(), ext_output.First(), double( 0));
          const size_t extxy2( ext_atom_map.First() * ext_atom_map.Second());
          size_t count( 0);
          for
          (
            const double *density_val( atom_map.Begin()), *density_val_end( atom_map.End());
            density_val != density_val_end;
            ++density_val, ++count
          )
          {
            if( *density_val != 0.0)
            {
              int indv = int( count);
              int k = indv / extxy2;
              indv -= k * extxy2;
              int j = indv / ext_atom_map.First();
              int i = indv - j * ext_atom_map.First();
              for( size_t indz( 0); indz < kernel_map.NumberLayers(); ++indz)
              {
                for( size_t indy( 0); indy < kernel_map.GetNumberRows(); ++indy)
                {
                  for( size_t indx( 0); indx < kernel_map.GetNumberCols(); ++indx)
                  {
                    output_map( k + indz, j + indy, i + indx) += kernel_map( indz, indy, indx) * ( *density_val);
                  }
                }
              }
            }
          }

//          if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Verbose))
//          {
//            std::ostringstream msg;
//            msg << "Lattice smoothing (sigma = atom rmsd) done: X Y Z "
//                << m_GridSpacing.X() * sqrt(varp) << ' ' << m_GridSpacing.Y() * sqrt(varp) << ' ' << m_GridSpacing.Z() * sqrt(varp)
//                << " Angstroem" << '\n'
//                << "Spatial resolution (2 sigma) of output map: ";
//            if( corrmode == 2)
//            {
//              const double reso( 2 * sqrt( math::Sqr( sigma) + varp * m_GridSpacing.X() * m_GridSpacing.X())); // variances are additive for uncorrelated samples
//              msg << reso << " slightly larger than target resolution due to uncorrected lattice smoothing";
//            }
//            else
//            {
//              msg << m_Resolution;
//            }
//            BCL_MessageVrb( msg.str());
//          }

          // index of map
          const linal::VectorND< int, 3> index
          (
            int( floor( grid3.X() / m_GridSpacing.X() + 0.5)),
            int( floor( grid3.Y() / m_GridSpacing.Y() + 0.5)),
            int( floor( grid3.Z() / m_GridSpacing.Z() + 0.5))
          );

          // intervals
          linal::VectorND< int, 3> intervals;
          for( size_t i( 0); i < 3; ++i)
          {
            if( index( i) <= -int( ext_output( i)))
            {
              intervals( i) = math::Absolute( index( i));
            }
            else if( index( i) >= 0)
            {
              intervals( i) = index( i) + ext_output( i) - 1;
            }
            else
            {
              intervals( i) = ext_output( i) - 1;
            }
          }

          // length
          const linal::Vector3D length
          (
            m_GridSpacing.X() * double( intervals( 0)),
            m_GridSpacing.Y() * double( intervals( 1)),
            m_GridSpacing.Z() * double( intervals( 2))
          );

          // construct map and return
          return Map
            (
              output_map,
              index,
              intervals,
              length,
              m_GridSpacing,
              Map::GetDefaultAngle(),
              Map::GetDefaultAxis(),
              linal::Vector3D( 0.0)
            );
        }

        default:
        {
          return Map();
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SimulateDefault::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_GridSpacing    , ISTREAM);
      io::Serialize::Read( m_Resolution     , ISTREAM);
      io::Serialize::Read( m_Margin         , ISTREAM);
      io::Serialize::Read( m_SmoothingKernel, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SimulateDefault::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_GridSpacing    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Resolution     , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Margin         , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SmoothingKernel, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace density
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
#include "density/bcl_density_simulate_gaussian_sphere.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SimulateGaussianSphere::s_Instance
    (
      GetObjectInstances().AddInstance( new SimulateGaussianSphere( linal::Vector3D( 2.3, 2.3, 2.3), 6.9))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from grid spacing and resolution
    //! @param GRID_SPACING the spacing for the density grid
    //! @param RESOLUTION the resolution to simulate for
    SimulateGaussianSphere::SimulateGaussianSphere( const linal::Vector3D &GRID_SPACING, const double RESOLUTION) :
      m_GridSpacing( GRID_SPACING),
      m_Resolution( RESOLUTION),
      m_Margin( 2)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SimulateGaussianSphere
    SimulateGaussianSphere *SimulateGaussianSphere::Clone() const
    {
      return new SimulateGaussianSphere( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SimulateGaussianSphere::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the resolution
    //! @param RESOLUTION the resolution for the density map to be generated
    void SimulateGaussianSphere::SetResolution( const double RESOLUTION)
    {
      m_Resolution = RESOLUTION;
    }

    //! @brief set the resolution
    double SimulateGaussianSphere::GetResolution() const
    {
      return m_Resolution;
    }

    //! @brief set the grid spacing
    //! @param GRID_SPACING the width of a grid element in x, y and z
    void SimulateGaussianSphere::SetGridSpacing( const linal::Vector3D &GRID_SPACING)
    {
      m_GridSpacing = GRID_SPACING;
    }

    //! @brief set the margin
    //! @param MARGIN number of additional cells next to last atom occupied cells
    void SimulateGaussianSphere::SetMargin( const size_t MARGIN)
    {
      m_Margin = MARGIN;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate simulated density from given list of atoms
    //! @param ATOMS siptrvector of atoms
    //! @return a simulated density map
    Map SimulateGaussianSphere::operator()( const util::SiPtrVector< const biol::Atom> &ATOMS) const
    {
      // measure ATOMS extent
      storage::VectorND< 2, linal::Vector3D> min_max_coord( DetermineGridCorners( ATOMS));

      linal::Vector3D &mincoord( min_max_coord.First());
      linal::Vector3D &maxcoord( min_max_coord.Second());

      // for undefined corners, return empty map
      if( !mincoord.IsDefined() || !maxcoord.IsDefined())
      {
        return Map();
      }

      // add margin
      mincoord -= m_Margin * m_GridSpacing;
      maxcoord += m_Margin * m_GridSpacing;

      // determine index
      const linal::VectorND< int, 3> index
      (
        int( std::floor( mincoord.X() / m_GridSpacing.X())),
        int( std::floor( mincoord.Y() / m_GridSpacing.Y())),
        int( std::floor( mincoord.Z() / m_GridSpacing.Z()))
      );

      // dimensions of grid
      const storage::VectorND< 3, size_t> dimensions
      (
        size_t( std::ceil( maxcoord.X() / m_GridSpacing.X())) - index( 0) ,
        size_t( std::ceil( maxcoord.Y() / m_GridSpacing.Y())) - index( 1),
        size_t( std::ceil( maxcoord.Z() / m_GridSpacing.Z())) - index( 2)
      );
      // allocate mask size and set values to 0
      math::Tensor< double> grid
      (
        dimensions.Third(), dimensions.Second(), dimensions.First(), double( 0.0)
      );

      // store the real space index for easier access
      const linal::Vector3D realspaceindex
        (
          index( 0)  * m_GridSpacing.X(),
          index( 1) * m_GridSpacing.Y(),
          index( 2)  * m_GridSpacing.Z()
        );

      // constants describing gaussian blob shape
      const double blob_k( math::Sqr( math::g_Pi / ( 2.4 + 0.8 * m_Resolution)));
      const double blob_c( math::Pow( blob_k / math::g_Pi, 1.5));

      // cutoff distance at 3 sigma
      const double cutoff_square( math::Sqr( 3.0 * ( 1.0 / math::Sqrt( 2.0)) * ( ( 2.4 + 0.8 * m_Resolution) / math::g_Pi)));
      const double cutoff( math::Sqrt( cutoff_square));

      // iterate over all atoms
      for
      (
        util::SiPtrVector< const biol::Atom>::const_iterator itr( ATOMS.Begin()), itr_end( ATOMS.End());
        itr != itr_end;
        ++itr
      )
      {
        const linal::Vector3D &current_coord( ( *itr)->GetCoordinates());

        // skip undefined coordinates
        if( !current_coord.IsDefined())
        {
          continue;
        }

        //  calculate the position within the grid
        const linal::Vector3D coord_in_grid( current_coord - realspaceindex);

        // mass of atom times c constant for sphere/blob
        const double blob_c_mass( blob_c * ( *itr)->GetType()->GetElementType()->GetProperty( chemistry::ElementTypeData::e_Mass));

        // min and max index within grid
        const storage::VectorND< 3, size_t> min_index
        (
          size_t( std::max( int( 0), int( std::floor( ( coord_in_grid.X() - cutoff) / m_GridSpacing.X())))),
          size_t( std::max( int( 0), int( std::floor( ( coord_in_grid.Y() - cutoff) / m_GridSpacing.Y())))),
          size_t( std::max( int( 0), int( std::floor( ( coord_in_grid.Z() - cutoff) / m_GridSpacing.Z()))))
        );
        const storage::VectorND< 3, size_t> max_index
        (
          size_t( std::min( int( dimensions.First()), int( std::ceil( ( coord_in_grid.X() + cutoff) / m_GridSpacing.X())))),
          size_t( std::min( int( dimensions.Second()), int( std::ceil( ( coord_in_grid.Y() + cutoff) / m_GridSpacing.Y())))),
          size_t( std::min( int( dimensions.Third()), int( std::ceil( ( coord_in_grid.Z() + cutoff) / m_GridSpacing.Z()))))
        );

        // iterate over voxel in grid to calculate the intensity/density
        for( size_t i( min_index.First()); i < max_index.First(); ++i)
        {
          linal::Vector3D pos_voxel;
          pos_voxel.X() = i * m_GridSpacing.X();
          for( size_t j( min_index.Second()); j < max_index.Second(); ++j)
          {
            pos_voxel.Y() = j * m_GridSpacing.Y();
            for( size_t k( min_index.Third()); k < max_index.Third(); ++k)
            {
              pos_voxel.Z() = k * m_GridSpacing.Z();

              const double square_distance( ( coord_in_grid - pos_voxel).SquareNorm());

              // check distance cutoff
              if( square_distance <= cutoff_square)
              {
                // calculate the distance between point in mask grid and specified coordinate
                grid( k, j, i) += blob_c_mass * std::exp( -blob_k * square_distance);
              }
            } // x
          } // y
        } // z
      } // iterate over atoms

      // intervals
      linal::VectorND< int, 3> intervals;
      for( size_t i( 0); i < 3; ++i)
      {
        if( index( i) <= -int( index( i)))
        {
          intervals( i) = math::Absolute( index( i));
        }
        else if( index( i) >= 0)
        {
          intervals( i) = index( i) + index( i) - 1;
        }
        else
        {
          intervals( i) = index( i) - 1;
        }
      }

      // length
      const linal::Vector3D length
      (
        m_GridSpacing.X() * double( intervals( 0)),
        m_GridSpacing.Y() * double( intervals( 1)),
        m_GridSpacing.Z() * double( intervals( 2))
      );

      // end
      return Map
             (
               grid,
               index,
               intervals,
               length,
               m_GridSpacing,
               Map::GetDefaultAngle(),
               Map::GetDefaultAxis(),
               linal::Vector3D( 0.0)
             );
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
    std::istream &SimulateGaussianSphere::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_GridSpacing    , ISTREAM);
      io::Serialize::Read( m_Resolution     , ISTREAM);
      io::Serialize::Read( m_Margin         , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SimulateGaussianSphere::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_GridSpacing    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Resolution     , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Margin         , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace density
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
#include "density/bcl_density_simulate_interface.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_mask_3d.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief measure the maximal map extend
    //! @param ATOMS list of atoms
    //! @return min coord xyz and max coord xyz
    storage::VectorND< 2, linal::Vector3D> SimulateInterface::DetermineGridCorners
    (
      const util::SiPtrVector< const biol::Atom> &ATOMS
    )
    {
      util::SiPtrVector< const linal::Vector3D> coordinates;
      coordinates.AllocateMemory( ATOMS.GetSize());

      // iterate over all atoms to acquire coordinates
      for
      (
        util::SiPtrVector< const biol::Atom>::const_iterator atom_itr( ATOMS.Begin()), atom_itr_end( ATOMS.End());
        atom_itr != atom_itr_end;
        ++atom_itr
      )
      {
        coordinates.PushBack( util::ToSiPtr( ( *atom_itr)->GetCoordinates()));
      }

      // determine grid corners for the coordinates
      return Mask3d::DetermineGridCorners( coordinates);
    }

  } // namespace density
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
#include "density/bcl_density_simulators.h"

// includes from bcl - sorted alphabetically
#include "density/bcl_density_simulate_default.h"
#include "density/bcl_density_simulate_gaussian_sphere.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

  //////////
  // data //
  //////////

    //! @brief default grid spacing
    const linal::Vector3D &Simulators::GetDefaultGridSpacing()
    {
      static const linal::Vector3D s_default_grid_spacing( 2.2, 2.2, 2.2);
      return s_default_grid_spacing;
    }

    //! @brief default resolution for density simulation
    double Simulators::GetDefaultResolution()
    {
      return 6.6;
    }

    //! @brief default constructor
    Simulators::Simulators() :
      e_GaussianSphere  ( AddEnum( "GaussianSphere"                        , util::ShPtr< SimulateInterface>( new SimulateGaussianSphere( GetDefaultGridSpacing(), GetDefaultResolution())))),
      e_NoSmoothing     ( AddEnum( "TrilinearInterpolation"                , util::ShPtr< SimulateInterface>( new SimulateDefault( GetDefaultGridSpacing(), GetDefaultResolution(), SimulateDefault::e_HardSphere)))),
      e_Gaussian        ( AddEnum( "TrilinearInterpolationGaussian"        , util::ShPtr< SimulateInterface>( new SimulateDefault( GetDefaultGridSpacing(), GetDefaultResolution(), SimulateDefault::e_Gaussian)))),
      e_Triangular      ( AddEnum( "TrilinearInterpolationTriangular"      , util::ShPtr< SimulateInterface>( new SimulateDefault( GetDefaultGridSpacing(), GetDefaultResolution(), SimulateDefault::e_Triangular)))),
      e_SemiEpanechnikov( AddEnum( "TrilinearInterpolationSemiEpanechnikov", util::ShPtr< SimulateInterface>( new SimulateDefault( GetDefaultGridSpacing(), GetDefaultResolution(), SimulateDefault::e_SemiEpanechnikov)))),
      e_Epanechnikov    ( AddEnum( "TrilinearInterpolationEpanechnikov"    , util::ShPtr< SimulateInterface>( new SimulateDefault( GetDefaultGridSpacing(), GetDefaultResolution(), SimulateDefault::e_Epanechnikov)))),
      e_HardSphere      ( AddEnum( "TrilinearInterpolationHardSphere"      , util::ShPtr< SimulateInterface>( new SimulateDefault( GetDefaultGridSpacing(), GetDefaultResolution(), SimulateDefault::e_HardSphere))))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Simulators::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief construct a simulator from given grid spacing and resolution
    //! @param SIMULATOR the simluator enum
    //! @param RESOLUTION desired resolution for density map
    //! @brief GRID_SPACING desired spacing for density map
    //! @return ShPtr to SimulatInterface for that simluator with grid spacing and resolution set
    util::ShPtr< SimulateInterface> Simulators::CreateSimulator
    (
      const Simulator &SIMULATOR,
      const linal::Vector3D &GRID_SPACING,
      const double RESOLUTION
    ) const
    {
      // copy the SimulateInterface derived class for that enumerator
      util::ShPtr< SimulateInterface> simulator( SIMULATOR->HardCopy());

      // set the grid spacing and resolution
      simulator->SetGridSpacing( GRID_SPACING);
      simulator->SetResolution( RESOLUTION);

      // end
      return simulator;
    }

    //! @brief construct on access function for all Simulators
    //! @return reference to only instances of Simulators
    Simulators &GetSimulators()
    {
      return Simulators::GetEnums();
    }

  } // namespace density

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< density::SimulateInterface>, density::Simulators>;

  } // namespace util
} // namespace bcl
