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

#ifndef BCL_COORD_POINT_CLOUD_H_
#define BCL_COORD_POINT_CLOUD_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_coord_movable_interface.h"
#include "linal/bcl_linal_vector_3d.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PointCloud
    //! @brief template class PointCloud
    //! @details TODO: document
    //!
    //! @see @link example_coord_point_cloud.cpp @endlink
    //! @author woetzen
    //! @date Oct 29, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PointCloud :
      public MovableInterface
    {

    private:

    //////////
    // data //
    //////////

      storage::Vector< linal::Vector3D> m_Data;

    public:

      typedef storage::Vector< linal::Vector3D>::iterator iterator;
      typedef storage::Vector< linal::Vector3D>::const_iterator const_iterator;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      PointCloud()
      {
      }

      //! construct from points
      PointCloud( const storage::Vector< linal::Vector3D> &POINTS) :
        m_Data( POINTS)
      {
      }

      //! copy constructor
      PointCloud *Clone() const
      {
        return new PointCloud( *this);
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

      size_t GetSize() const
      {
        return m_Data.GetSize();
      }

      iterator Begin()
      {
        return m_Data.Begin();
      }

      iterator End()
      {
        return m_Data.End();
      }

      const_iterator Begin() const
      {
        return m_Data.Begin();
      }

      const_iterator End() const
      {
        return m_Data.End();
      }

      void PushBack( const linal::Vector3D &POINT)
      {
        m_Data.PushBack( POINT);
      }

      operator const storage::Vector< linal::Vector3D> &() const
      {
        return m_Data;
      }

      operator storage::Vector< linal::Vector3D> &()
      {
        return m_Data;
      }

      const storage::Vector< linal::Vector3D> &GetData() const
      {
        return m_Data;
      }

      storage::Vector< linal::Vector3D> &GetData()
      {
        return m_Data;
      }

    ////////////////
    // operations //
    ////////////////

      //! returns the geometric center of the object
      linal::Vector3D GetCenter() const
      {
        return CenterOfMass( util::ConvertToConstSiPtrVector( m_Data));
      }

      //! translate the object along a given TRANSLATION vector
      void Translate( const linal::Vector3D &TRANSLATION)
      {
        for( iterator point_itr( Begin()), point_itr_end( End()); point_itr != point_itr_end; ++point_itr)
        {
          point_itr->Translate( TRANSLATION);
        }
      }

      //! transform the object by a given TRANSFORMATIONMATRIX3D
      void Transform( const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D)
      {
        for( iterator point_itr( Begin()), point_itr_end( End()); point_itr != point_itr_end; ++point_itr)
        {
          point_itr->Transform( TRANSFORMATIONMATRIX3D);
        }
      }

      //! rotate the object by a given ROTATIONMATRIX3D
      void Rotate( const math::RotationMatrix3D &ROTATIONMATRIX3D)
      {
        for( iterator point_itr( Begin()), point_itr_end( End()); point_itr != point_itr_end; ++point_itr)
        {
          point_itr->Rotate( ROTATIONMATRIX3D);
        }
      }

    ////////////////
    // operations //
    ////////////////

      //! Translate the entire pointcloud to the center of mass
      void CenterPointCloud();

      //! @brief returns a mpa that maps all points, to al list of coordinates, that are neighbors of that point within the given distance
      storage::Map< util::SiPtr< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
      DetermineNeighbors( const double DISTANCE) const;

      //! @brief remove all points from the cloud that do not have a certain number of neighbors within a distance
      //! @param MIN_NEIGHBORS minimal number of neighbors
      //! @param DISTANCE distance in which the minimal number of neighbors have to be found
      size_t RemoveSingles( const size_t MIN_NEIGHBORS, const double DISTANCE);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! write this PointCloud to a given FILE
      //! @param FILENAME
      //! pdb file where pointcloud should be written to
      void WriteToPDB( const std::string &FILENAME) const;

    }; // class PointCloud

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_POINT_CLOUD_H_
