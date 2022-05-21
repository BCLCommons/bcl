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

#ifndef BCL_FIND_LOCATOR_COORDINATES_TRIGONAL_H_
#define BCL_FIND_LOCATOR_COORDINATES_TRIGONAL_H_

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
    //! @class LocatorCoordinatesTrigonal
    //! @brief for locating the coordinates of a point on the vertex of trigonal geometry
    //! @details The location of the coordinates is defined by two other vertices, a center point, and a distance from
    //! the center point. The calculated coordinates are perpendicular to the initial two vertices and does not
    //! necessarily exist in t_ArgumentType.
    //!
    //! @see @link example_find_locator_coordinates_trigonal.cpp @endlink
    //! @author alexanns
    //! @date Apr 12, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType>
    class LocatorCoordinatesTrigonal :
      public LocatorInterface< linal::Vector3D, t_ArgumentType>
    {

    private:

    //////////
    // data //
    //////////

      //! "m_Center" is a ShPtr to a LocatorInterface which takes t_ArgumentType and returns the
      //! MovableInterface whose coordinates will be used as the center point in defining the Trigonal geometry
      util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > m_Center;

      //! "m_VertexA" is a ShPtr to a LocatorInterface which takes t_ArgumentType and returns the MovableInterface
      //! whose coordinates will be used as one of the two vertices defining the Trigonal geometry
      util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > m_VertexA;

      //! "m_VertexB" is a ShPtr to a LocatorInterface which takes t_ArgumentType and returns the MovableInterface
      //! whose coordinates will be used as one of the two vertices defining the Trigonal geometry
      util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > m_VertexB;

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
      LocatorCoordinatesTrigonal();

      //! @brief constructor taking all member objects
      //! @param CENTER ShPtr< LocatorInterface> giving MoveableInterface with coords at center of trigonal geometry
      //! @param VERTEX_A ShPtr< LocatorInterface> giving MoveableInterface with coords at vertex of trigonal geometry
      //! @param VERTEX_B ShPtr< LocatorInterface> giving MoveableInterface with coords at vertex of trigonal geometry
      //! @param DISTANCE double denoting how far away from the center the located coordinates will be
      LocatorCoordinatesTrigonal
      (
        const util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > &CENTER,
        const util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > &VERTEX_A,
        const util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > &VERTEX_B,
        const double DISTANCE
      );

      //! @brief Clone function
      //! @return pointer to new LocatorCoordinatesTrigonal
      LocatorCoordinatesTrigonal< t_ArgumentType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Locate gives the coordinates of the the point defined by trigonal geometry according to the members
      //! @param ARGUMENT t_ArgumentType from which the coordinates will be located
      //! @return the coordinates of the point determined - does not necessarily exist in ARGUMENT
      linal::Vector3D Locate( const t_ArgumentType &ARGUMENT) const;

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
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // template class LocatorCoordinatesTrigonal

  //////////
  // data //
  //////////

    //! single instance of that class
    template< typename t_ArgumentType>
    const util::SiPtr< const util::ObjectInterface> LocatorCoordinatesTrigonal< t_ArgumentType>::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorCoordinatesTrigonal< t_ArgumentType>())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_ArgumentType>
    LocatorCoordinatesTrigonal< t_ArgumentType>::LocatorCoordinatesTrigonal() :
      m_Center(),
      m_VertexA(),
      m_VertexB(),
      m_DistanceFromCenter()
    {
    }

    //! @brief constructor taking all member objects
    //! @param CENTER ShPtr< LocatorInterface> giving MoveableInterface with coords at center of trigonal geometry
    //! @param VERTEX_A ShPtr< LocatorInterface> giving MoveableInterface with coords at vertex of trigonal geometry
    //! @param VERTEX_B ShPtr< LocatorInterface> giving MoveableInterface with coords at vertex of trigonal geometry
    //! @param DISTANCE double denoting how far away from the center the located coordinates will be
    template< typename t_ArgumentType>
    LocatorCoordinatesTrigonal< t_ArgumentType>::LocatorCoordinatesTrigonal
    (
      const util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > &CENTER,
      const util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > &VERTEX_A,
      const util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > &VERTEX_B,
      const double DISTANCE
    ) :
      m_Center( CENTER),
      m_VertexA( VERTEX_A),
      m_VertexB( VERTEX_B),
      m_DistanceFromCenter( DISTANCE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LocatorCoordinatesTrigonal
    template< typename t_ArgumentType>
    LocatorCoordinatesTrigonal< t_ArgumentType> *LocatorCoordinatesTrigonal< t_ArgumentType>::Clone() const
    {
      return new LocatorCoordinatesTrigonal( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_ArgumentType>
    const std::string &LocatorCoordinatesTrigonal< t_ArgumentType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Locate gives the coordinates of the the point defined by trigonal geometry according to the members
    //! @param ARGUMENT t_ArgumentType from which the coordinates will be located
    //! @return the coordinates of the point determined - does not necessarily exist in ARGUMENT
    template< typename t_ArgumentType> linal::Vector3D
    LocatorCoordinatesTrigonal< t_ArgumentType>::Locate( const t_ArgumentType &ARGUMENT) const
    {
      // locate the moveable interfaces that defined the vertices
      const linal::Vector3D center( m_Center->Locate( ARGUMENT));
      const linal::Vector3D vertex_a( m_VertexA->Locate( ARGUMENT));
      const linal::Vector3D vertex_b( m_VertexB->Locate( ARGUMENT));

      // true if any of the located movable interfaces SiPtrs is undefined is undefined indicating the movable interface
      // could not be found
      if( !center.IsDefined() || !vertex_a.IsDefined() || !vertex_b.IsDefined())
      {
        // return an undefined Vector3D indicating that location of the MoveableInterface failed
        return linal::Vector3D( util::GetUndefined< double>());
      }

      // get the coordinates of the point defining the trigonal geometry
      const linal::Vector3D trigonal_coordinates
      (
        linal::CoordinatesTrigonal( center, vertex_a, vertex_b, m_DistanceFromCenter)
      );

      // return the coordinates of the point defining the trigonal geometry
      return trigonal_coordinates;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_ArgumentType>
    std::istream &LocatorCoordinatesTrigonal< t_ArgumentType>::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Center            , ISTREAM);
      io::Serialize::Read( m_VertexA           , ISTREAM);
      io::Serialize::Read( m_VertexB           , ISTREAM);
      io::Serialize::Read( m_DistanceFromCenter, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    template< typename t_ArgumentType>
    std::ostream &LocatorCoordinatesTrigonal< t_ArgumentType>::Write
    (
      std::ostream &OSTREAM, const size_t INDENT
    ) const
    {
      // write member
      io::Serialize::Write( m_Center            , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_VertexA           , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_VertexB           , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DistanceFromCenter, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace find
} // namespace bcl

#endif // BCL_FIND_LOCATOR_COORDINATES_TRIGONAL_H_ 
