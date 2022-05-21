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

#ifndef BCL_DENSITY_CONNECTIVITY_H_
#define BCL_DENSITY_CONNECTIVITY_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_density_map.h"
#include "assemble/bcl_assemble_sse_geometry.h"
#include "coord/bcl_coord_geometry_interface.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Connectivity
    //! @brief has a collection of functions that will handle the calculation for connectivities of points
    //!        within density maps
    //!
    //! @see @link example_density_connectivity.cpp @endlink
    //! @author woetzen, linders
    //! @date 09/13/2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Connectivity :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      util::ShPtr< assemble::SSEGeometryInterface> m_BodyA;        //!< body of first density rod
      bool                                   m_BodyABegin;   //!< begin or end of body A used as point for connectivity
                                                             //!< (true for begin)
      util::ShPtr< assemble::SSEGeometryInterface> m_BodyB;        //!< body of second density rod
      bool                                   m_BodyBBegin;   //!< begin or end of body B used as point for connectivity
                                                             //!< (true for begin)
      double                                 m_Connectivity; //!< intensity at which points are connected
      double                                 m_Distance;     //!< real space distance between two points

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Connectivity();

      //! @brief Connectivity from various arguments
      Connectivity
      (
        const util::ShPtr< assemble::SSEGeometryInterface> &BODY_A,
        const bool BODY_A_BEGIN,
        const util::ShPtr< assemble::SSEGeometryInterface> &BODY_B,
        const bool BODY_B_BEGIN,
        const double CONNECTIVITY = util::GetUndefinedDouble(),
        const double DISTANCE = util::GetUndefinedDouble()
      );

      //!virtual copy constructor
      Connectivity *Clone() const
      {
        return new Connectivity( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return connectivity of Connectivity object
      //! @return the connectivity of Connectivity object as double
      double GetConnectivity() const
      {
        return m_Connectivity;
      }

      //! @brief return distance of Connectivity object
      //! @return the distance of Connectivity object as double
      double GetDistance() const
      {
        return m_Distance;
      }

      //! @brief return first body of Connectivity object
      //! @return a ShPtr to the first body of Connectivity object
      const util::ShPtr< assemble::SSEGeometryInterface> &GetBodyA() const
      {
        return m_BodyA;
      }

      //! @brief return second body of Connectivity object
      //! @return a ShPtr to the second body of Connectivity object
      const util::ShPtr< assemble::SSEGeometryInterface> &GetBodyB() const
      {
        return m_BodyB;
      }

      //! @brief return orientation of first body of Connectivity object
      //! @return the orientation of first body of Connectivity object as bool (true if begin, false if end)
      bool GetOrientationA() const
      {
        return m_BodyABegin;
      }

      //! @brief return orientation of second body of Connectivity object
      //! @return the orientation of second body of Connectivity object as bool (true if begin, false if end)
      bool GetOrientationB() const
      {
        return m_BodyBBegin;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief determine all 4 DensityConnectivities from every set of two bodies in a list in a density map
      //! @param DENSITY_MAP the density map the bodies are lying in
      //! @param BODIES list of bodies for which the mutual connectivities are to be determined
      //! @return list of connectivities of combined begin-begin, begin-end, end-begin, end-end for each pair of bodies
      static
      storage::List< Connectivity>
      DetermineConnectivities
      (
        const Map &DENSITY_MAP,
        const util::ShPtrList< assemble::SSEGeometryInterface> &BODIES
      );

      //! @brief determine all 4 DensityConnectivities from 2 bodies in a density map given an optional length reduction
      //! @param DENSITY_MAP the density map the bodies are lying in
      //! @param BODY_A body a
      //! @param BODY_B body b
      //! @return list of connectivities of combined begin-begin, begin-end, end-begin, end-end of both bodies
      static
      storage::List< Connectivity>
      DetermineConnectivities
      (
        const Map &DENSITY_MAP,
        const util::ShPtr< assemble::SSEGeometryInterface> &BODY_A,
        const util::ShPtr< assemble::SSEGeometryInterface> &BODY_B
      );

      //! @brief determine the lowest intensity the is necessary to connect two points in the electron density map
      //!        (this is the weakest point in the connection between two points)
      //! @param DENSITY_MAP the density map where the points are in
      //! @param COORDINATE_A the first point
      //! @param COORDINATE_B the second point
      //! @return an intensity that is low enough to connect the two points through voxels of intensities larger than the returned intensity
      static
      double
      ConnectionIntensity
      (
        const Map &DENSITY_MAP,
        const linal::Vector3D &COORDINATE_A,
        const linal::Vector3D &COORDINATE_B
      );

    private:

      //! @brief creates a subdensity that encapsulates the two coordinates in a 2 times the distance between the points large densitymap
      //! @param DENSITY_MAP the density map where the points are in
      //! @param COORDINATE_A the first point
      //! @param COORDINATE_B the second point
      //! @return the subdensity which contains the two points with generous spacing
      static
      Map
      SubDensityFromTwoCoordinates
      (
        const Map &DENSITY_MAP,
        const linal::Vector3D &COORDINATE_A,
        const linal::Vector3D &COORDINATE_B
      );

      //! @brief set the voxel for the current point to "nan" if its intensity is larger than the threshold
      //! @param DENSITY_MAP the map that contains the voxel
      //! @param THRESHOLD threshold for the intensity that determines if the point should be colored or not
      //! @param CURRENT_POINT the index of the voxel to be checked
      //! @return false, if it was already colored, if it is smaller than the threshold. True if it was colored.
      static
      bool
      ColorVoxel
      (
        Map &DENSITY_MAP,
        const double THRESHOLD,
        const linal::Vector< int> &CURRENT_POINT
      );

      //! @brief iterate over all neighboring voxels of the MIDDLE_POINT
      //! @param DENSITY_MAP map where the MIDDLE_POINT lies in
      //! @param THRESHOLD threshold for the coloring (all voxels that have intensity THRESHOLD and higher are colored)
      //! @param MIDDLE_POINT the center point from which the surrounding voxels are explored and colored if above THRESHOLD
      static
      void
      ColorNeighboringVoxels
      (
        Map &DENSITY_MAP,
        const double THRESHOLD,
        const linal::Vector< int> &MIDDLE_POINT
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read Connectivity from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! write Connectivity into std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class LessThan
      //! @brief orders density::Connectivity objects by sh-ptr addresses of the bodies, then by orientation
      //!
      //! @see @link example_density_connectivity.cpp @endlink
      //! @author woetzen, linders
      //! @date 09/13/2007
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct LessThan
      {
        //! @brief perform smaller than  comparison of Connectivity objects
        //! @param CONNECTIVITY_LHS  left hand argument
        //! @param CONNECTIVITY_RHS right hand argument
        //! @return bool - compares the two bodies on LHS and RHS, if they are equal, the bools describing the beginning
        //!         of A and B are compared
        bool operator()
        (
          const Connectivity &CONNECTIVITY_LHS,
          const Connectivity &CONNECTIVITY_RHS
        ) const
        {
          if( CONNECTIVITY_LHS.GetBodyA() == CONNECTIVITY_RHS.GetBodyA())
          {
            if( CONNECTIVITY_LHS.GetBodyB() == CONNECTIVITY_RHS.GetBodyB())
            {
              if
              (
                CONNECTIVITY_LHS.GetOrientationA() == CONNECTIVITY_RHS.GetOrientationA()
              )
              {
                if
                (
                  CONNECTIVITY_LHS.GetOrientationB() == CONNECTIVITY_RHS.GetOrientationB()
                )
                {
                  return false;
                }
                return CONNECTIVITY_LHS.GetOrientationB() < CONNECTIVITY_RHS.GetOrientationB();
              }
              return CONNECTIVITY_LHS.GetOrientationA() < CONNECTIVITY_RHS.GetOrientationA();
            }
            return CONNECTIVITY_LHS.GetBodyB() < CONNECTIVITY_RHS.GetBodyB();
          }
          return CONNECTIVITY_LHS.GetBodyA() < CONNECTIVITY_RHS.GetBodyA();
        }
      }; // struct LessThan

    }; // class Connectivity

    //! @brief boolean operator CONNECTIVITY_LHS == CONNECTIVITY_RHS
    //! @param CONNECTIVITY_LHS first CONNECTIVITY
    //! @param CONNECTIVITY_RHS second CONNECTIVITY
    //! @return whether CONNECTIVITY_LHS is equal to CONNECTIVITY_RHS
    BCL_API bool operator ==( const Connectivity &CONNECTIVITY_LHS, const Connectivity &CONNECTIVITY_RHS);

  } // namespace density
} // namespace bcl

#endif //BCL_DENSITY_CONNECTIVITY_H_
