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

#ifndef BCL_FIND_LOCATOR_COORDINATES_TETRAHEDRAL_H_
#define BCL_FIND_LOCATOR_COORDINATES_TETRAHEDRAL_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_find_locator_interface.h"
#include "coord/bcl_coord_movable_interface.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorCoordinatesTetrahedral
    //! @brief for locating the coordinates of a point on the vertex of a tetrahedron
    //! @details The location of the coordinates is defined by three other vertices, a center point, and a distance from
    //! the center point So the located coordinates are not strictly defined by a tetrahedron because the distances
    //! between the three defined vertices won't necessarily be the same distance as from them to the located
    //! coordinates
    //!
    //! @see @link example_find_locator_coordinates_tetrahedral.cpp @endlink
    //! @author alexanns
    //! @date Apr 12, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType>
    class LocatorCoordinatesTetrahedral :
      public LocatorInterface< linal::Vector3D, t_ArgumentType>
    {

    private:

    //////////
    // data //
    //////////

      //! "m_Center" is a ShPtr to a LocatorInterface which takes t_ArgumentType and returns the
      //! MovableInterface whose coordinates will be used as the center point in defining the tetrahedral geometry
      util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > m_Center;

      //! "m_VertexA" is a ShPtr to a LocatorInterface which takes t_ArgumentType and returns the MovableInterface
      //! whose coordinates will be used as one of the three vertices defining the tetrahedral geometry
      util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > m_VertexA;

      //! "m_VertexB" is a ShPtr to a LocatorInterface which takes t_ArgumentType and returns the MovableInterface
      //! whose coordinates will be used as one of the three vertices defining the tetrahedral geometry
      util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > m_VertexB;

      //! "m_VertexC" is a ShPtr to a LocatorInterface which takes t_ArgumentType and returns the MovableInterface
      //! whose coordinates will be used as one of the three vertices defining the tetrahedral geometry
      util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > m_VertexC;

      //! "m_DistanceFromCenter" is a double denoting how far away from the center the located coordinates will be
      double m_DistanceFromCenter;

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
      LocatorCoordinatesTetrahedral() :
        m_Center(),
        m_VertexA(),
        m_VertexB(),
        m_VertexC(),
        m_DistanceFromCenter()
      {
      }

      //! @brief constructor taking all member objects
      //! @param CENTER ShPtr< LocatorInterface> giving MoveableInterface with coords at center of tetrahedron
      //! @param VERTEX_A ShPtr< LocatorInterface> giving MoveableInterface with coords at vertex of tetrahedron
      //! @param VERTEX_B ShPtr< LocatorInterface> giving MoveableInterface with coords at vertex of tetrahedron
      //! @param VERTEX_C ShPtr< LocatorInterface> giving MoveableInterface with coords at vertex of tetrahedron
      //! @param DISTANCE double denoting how far away from the center the located coordinates will be
      LocatorCoordinatesTetrahedral
      (
        const LocatorInterface< linal::Vector3D, t_ArgumentType> &CENTER,
        const LocatorInterface< linal::Vector3D, t_ArgumentType> &VERTEX_A,
        const LocatorInterface< linal::Vector3D, t_ArgumentType> &VERTEX_B,
        const LocatorInterface< linal::Vector3D, t_ArgumentType> &VERTEX_C,
        const double DISTANCE
      ) :
        m_Center( CENTER.Clone()),
        m_VertexA( VERTEX_A.Clone()),
        m_VertexB( VERTEX_B.Clone()),
        m_VertexC( VERTEX_C.Clone()),
        m_DistanceFromCenter( DISTANCE)
      {
      }

      //! @brief constructor taking all member objects
      //! @param CENTER ShPtr< LocatorInterface> giving MoveableInterface with coords at center of tetrahedron
      //! @param VERTEX_A ShPtr< LocatorInterface> giving MoveableInterface with coords at vertex of tetrahedron
      //! @param VERTEX_B ShPtr< LocatorInterface> giving MoveableInterface with coords at vertex of tetrahedron
      //! @param VERTEX_C ShPtr< LocatorInterface> giving MoveableInterface with coords at vertex of tetrahedron
      //! @param DISTANCE double denoting how far away from the center the located coordinates will be
      LocatorCoordinatesTetrahedral
      (
        const util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > &CENTER,
        const util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > &VERTEX_A,
        const util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > &VERTEX_B,
        const util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > &VERTEX_C,
        const double DISTANCE
      ) :
        m_Center( CENTER),
        m_VertexA( VERTEX_A),
        m_VertexB( VERTEX_B),
        m_VertexC( VERTEX_C),
        m_DistanceFromCenter( DISTANCE)
      {
      }

      //! @brief Clone function
      //! @return pointer to new LocatorCoordinatesTetrahedral
      LocatorCoordinatesTetrahedral< t_ArgumentType> *Clone() const
      {
        return new LocatorCoordinatesTetrahedral( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Locate gives the coordinates of the the point defined by tetrahedral geometry according to the members
      //! @param ARGUMENT t_ArgumentType from which the coordinates will be located
      //! @return the coordinates of the point determined
      linal::Vector3D Locate( const t_ArgumentType &ARGUMENT) const
      {
        // locate the moveable interfaces that defined the vertices
        const linal::Vector3D center( m_Center->Locate( ARGUMENT));
        const linal::Vector3D vertex_a( m_VertexA->Locate( ARGUMENT));
        const linal::Vector3D vertex_b( m_VertexB->Locate( ARGUMENT));
        const linal::Vector3D vertex_c( m_VertexC->Locate( ARGUMENT));

        // true if any of the located movable interfaces SiPtrs is undefined
        // indicating the movable interface could not be found
        if( !center.IsDefined() || !vertex_a.IsDefined() || !vertex_b.IsDefined() || !vertex_c.IsDefined())
        {
          // return an undefined Vector3D indicating that location of the MoveableInterface failed
          return linal::Vector3D( util::GetUndefined< double>());
        }

        return linal::CoordinatesTetrahedral( center, vertex_a, vertex_b, vertex_c, m_DistanceFromCenter);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read member
        io::Serialize::Read( m_Center, ISTREAM);
        io::Serialize::Read( m_VertexA, ISTREAM);
        io::Serialize::Read( m_VertexB, ISTREAM);
        io::Serialize::Read( m_VertexC, ISTREAM);
        io::Serialize::Read( m_DistanceFromCenter, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_Center, OSTREAM, INDENT)  << '\n';
        io::Serialize::Write( m_VertexA, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_VertexB, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_VertexC, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_DistanceFromCenter, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class LocatorCoordinatesTetrahedral

    //! single instance of that class
    template< typename t_ArgumentType>
    const util::SiPtr< const util::ObjectInterface> LocatorCoordinatesTetrahedral< t_ArgumentType>::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorCoordinatesTetrahedral< t_ArgumentType>())
    );

  } // namespace find
} // namespace bcl

#endif // BCL_FIND_LOCATOR_COORDINATES_TETRAHEDRAL_H_ 
