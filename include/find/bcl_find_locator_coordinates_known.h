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

#ifndef BCL_FIND_LOCATOR_COORDINATES_KNOWN_H_
#define BCL_FIND_LOCATOR_COORDINATES_KNOWN_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_find_locator_interface.h"
#include "coord/bcl_coord_movable_interface.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorCoordinatesKnown
    //! @brief Class for locating the center coordinates of a movable interface object that is specified by a locator.
    //! @details This class uses a locator to locate a movable interface object then returns the coordinates of the
    //!          center of the movable interface object. This class is important for code where a LocatorInterface
    //!          is needed that returns coordinates especially in classes related to nmr restraints. This class is
    //!          related to the LocatorCoordinatesTrigonal and LocatorCoordinatesTetrahedral classes.
    //!
    //! @see @link example_find_locator_coordinates_known.cpp @endlink
    //! @author alexanns
    //! @date May 28, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType>
    class LocatorCoordinatesKnown :
      public LocatorInterface< linal::Vector3D, t_ArgumentType>
    {

    private:

    //////////
    // data //
    //////////

      //! "m_ObjectWithCoordinatesLocator" is a ShPtr to a LocatorInterface which takes t_ArgumentType and returns the
      //! MovableInterface whose coordinates are desired
      util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > m_ObjectWithCoordinatesLocator;

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
      LocatorCoordinatesKnown() :
        m_ObjectWithCoordinatesLocator()
      {
      }

      //! @brief constructor taking ShPtr to LocatorInterface which takes t_ArgumentType and returns a MovableInterface
      //! @param LOCATOR ShPtr to LocatorInterface which will give the MoveableInterface whose coordinates are desired
      LocatorCoordinatesKnown
      (
        const LocatorInterface< linal::Vector3D, t_ArgumentType> &LOCATOR
      ) :
        m_ObjectWithCoordinatesLocator( LOCATOR.Clone())
      {
      }

      //! @brief constructor taking ShPtr to LocatorInterface which takes t_ArgumentType and returns a MovableInterface
      //! @param LOCATOR ShPtr to LocatorInterface which will give the MoveableInterface whose coordinates are desired
      LocatorCoordinatesKnown
      (
        const util::ShPtr< LocatorInterface< linal::Vector3D, t_ArgumentType> > &LOCATOR
      ) :
        m_ObjectWithCoordinatesLocator( LOCATOR)
      {
      }

      //! @brief Clone function
      //! @return pointer to new LocatorCoordinatesKnown
      LocatorCoordinatesKnown< t_ArgumentType> *Clone() const
      {
        return new LocatorCoordinatesKnown( *this);
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

      //! @brief Locate gives the center coordinates of the MovableInterface whose coordinates are desired
      //! @param ARGUMENT t_ArgumentType from which the coordinates will be located
      //! @return the center coordinates of the MovableInterface located by "m_ObjectWithCoordinatesLocator"
      linal::Vector3D Locate( const t_ArgumentType &ARGUMENT) const
      {
        // create SiPtr to MovableInterface "located" initialize with what is found by "m_ObjectWithCoordinatesLocator"
        const linal::Vector3D located( m_ObjectWithCoordinatesLocator->Locate( ARGUMENT));

        // true if the SiPtr to MovableInterface is not defined indicating that it could not be located
        if( !located.IsDefined())
        {
          // return empty Vector3D to signify location was not successfull
          return linal::Vector3D();
        }

        // return the center coordinates of the MovableInterface that was located
        return located;
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
        io::Serialize::Read( m_ObjectWithCoordinatesLocator, ISTREAM);

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
        io::Serialize::Write( m_ObjectWithCoordinatesLocator, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class LocatorCoordinatesKnown

    //! single instance of that class
    template< typename t_ArgumentType>
    const util::SiPtr< const util::ObjectInterface> LocatorCoordinatesKnown< t_ArgumentType>::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorCoordinatesKnown< t_ArgumentType>())
    );

  } // namespace find
} // namespace bcl

#endif // BCL_FIND_LOCATOR_COORDINATES_KNOWN_H_ 
