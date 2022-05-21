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

#ifndef BCL_FIND_LOCATOR_COORDINATES_AVERAGE_H_
#define BCL_FIND_LOCATOR_COORDINATES_AVERAGE_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_find_locator_coordinates_interface.h"
#include "bcl_find_locator_coordinates_known.h"
#include "bcl_find_locator_interface.h"
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d.h"
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorCoordinatesAverage
    //! @brief Locates the average position of a list of coordinates
    //! @details Locates the average position of a list of located coordinates
    //!
    //! @tparam t_Argument the type of object from which the average coordinate will be located
    //!
    //! @see @link example_find_locator_coordinates_average.cpp @endlink
    //! @author alexanns
    //! @date Feb 14, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Argument>
    class LocatorCoordinatesAverage :
      public LocatorCoordinatesInterface< t_Argument>
    {

    private:

    //////////
    // data //
    //////////

      //! list of locators to find coordinates that will be averaged
      storage::List< util::Implementation< LocatorCoordinatesInterface< t_Argument> > > m_CoordinateLocators;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorCoordinatesAverage() :
        m_CoordinateLocators()
      {
        m_CoordinateLocators.PushBack
        (
          util::Implementation< LocatorCoordinatesInterface< t_Argument> >( assemble::LocatorAtom( 'A', util::GetUndefinedSize_t(), biol::GetAtomTypes().CA))
        );
        m_CoordinateLocators.PushBack
        (
          util::Implementation< LocatorCoordinatesInterface< t_Argument> >( assemble::LocatorAtom( 'A', util::GetUndefinedSize_t(), biol::GetAtomTypes().CB))
        );
      }

      //! @brief constructor taking member variable as parameter
      //! @param LOCATORS locators that will find coordinates to be averaged
      LocatorCoordinatesAverage( const storage::List< util::Implementation< LocatorCoordinatesInterface< t_Argument> > > &LOCATORS) :
        m_CoordinateLocators( LOCATORS)
      {
      }

      //! @brief Clone function
      //! @return pointer to new LocatorCoordinatesAverage< t_Argument>
      LocatorCoordinatesAverage< t_Argument> *Clone() const
      {
        return new LocatorCoordinatesAverage< t_Argument>( *this);
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

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const
      {
        static const std::string s_Name( "LocatorCoordinatesAverage");
        return s_Name;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief locates the average of the desired coordinates in ARGUMENT
      //! @param ARGUMENT the object where coordinates will be located and averaged
      //! @return linal::Vector3D which is the average of the coordinates located in ARGUMENT
      linal::Vector3D Locate( const t_Argument &ARGUMENT) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer parameters;
        parameters.SetClassDescription( "Locates the average position of a list of located coordinates.");
        parameters.AddInitializer
        (
          "",
          "methods for locating coordinates that will be averaged",
          io::Serialization::GetAgent( &m_CoordinateLocators)
        );

        return parameters;
      }

    }; // template class LocatorCoordinatesAverage

    // instantiate s_Instance
    template< typename t_Argument>
    const util::SiPtr< const util::ObjectInterface> LocatorCoordinatesAverage< t_Argument>::s_Instance
    (
      util::Enumerated< LocatorCoordinatesInterface< t_Argument> >::AddInstance( new LocatorCoordinatesAverage< t_Argument>())
    );

  ////////////////
  // operations //
  ////////////////

    //! @brief locates the average of the desired coordinates in ARGUMENT
    //! @param ARGUMENT the object where coordinates will be located and averaged
    //! @return linal::Vector3D which is the average of the coordinates located in ARGUMENT
    template< typename t_Argument>
    linal::Vector3D LocatorCoordinatesAverage< t_Argument>::Locate( const t_Argument &ARGUMENT) const
    {
      // will hold the average of the coordinats located in ARGUMENT
      linal::Vector3D average_coord( 0, 0, 0);

      // will hold the number of coordinates that are summed
      size_t num_coords( 0);

      // iterate through the locators in m_CoordinateLocators to sum up the coordinates
      for
      (
        typename storage::List< util::Implementation< LocatorCoordinatesInterface< t_Argument> > >::const_iterator
          itr( m_CoordinateLocators.Begin()), itr_end( m_CoordinateLocators.End());
        itr != itr_end; ++itr
      )
      {
        // locate the coordinates in ARGUMENT
        const linal::Vector3D current_coords( ( *itr)->Locate( ARGUMENT));

        // true if the coordinates are defined
        if( current_coords.IsDefined())
        {
          // add the coordinates to average_coord and increment num_coords
          average_coord += current_coords;
          ++num_coords;
        }
      }

      // true if no coordinates were located
      if( !num_coords)
      {
        // return undefined vector3D
        BCL_MessageDbg( "no coordinates located");
        return linal::Vector3D( util::GetUndefinedDouble(), util::GetUndefinedDouble(), util::GetUndefinedDouble());
      }

      // divide average coord by the number of coordinates that were added into it
      average_coord /= double( num_coords);

      // return the average coordinate of all that were located
      return average_coord;
    }

  } // namespace find

} // namespace bcl

#endif // BCL_FIND_LOCATOR_COORDINATES_AVERAGE_H_
