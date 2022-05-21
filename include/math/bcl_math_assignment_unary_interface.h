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

#ifndef BCL_MATH_ASSIGNMENT_UNARY_INTERFACE_H_
#define BCL_MATH_ASSIGNMENT_UNARY_INTERFACE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AssignmentUnaryInterface
    //! @details Interface for function objects that perform unary assignment operations ( Not(), Abs(), Log(), etc.)
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Mar 11, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AssignmentUnaryInterface :
      public util::SerializableInterface
    {

    public:

      //! @brief get all the different unary assignments
      //! @return vector of all different unary assignment types
      static const storage::Vector< util::Implementation< AssignmentUnaryInterface> > &GetUnaryAssignments();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      virtual AssignmentUnaryInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief operator the implements the operation on the argument
      //! @param VALUE the value to operate on
      virtual void operator()( float &VALUE) const = 0;

      //! @brief operator the implements the operation on the argument
      //! @param VALUE the value to operate on
      virtual void operator()( double &VALUE) const = 0;

      //! @brief helper function to perform the given assignment operation on all the items in a container
      //! @param CONTAINER the container of interest
      template< typename t_Container>
      void PerformOnEach( t_Container &CONTAINER)
      {
        for( typename t_Container::iterator itr( CONTAINER.Begin()), itr_end( CONTAINER.End()); itr != itr_end; ++itr)
        {
          operator()( *itr);
        }
      }

    public:

      //! @brief function to add all instances to the enumerated interface for
      //! @tparam t_InterfaceType the enumerated type to add all newly-created t_ClassTypes to
      //! @tparam t_ClassType the actual class instance; must take a single assignment in its constructor
      template< typename t_InterfaceType, typename t_ClassType>
      static util::ObjectInterface *AddUnaryInstances()
      {
        // keep a pointer to the last created instance
        util::ObjectInterface *last_instance( NULL);
        for
        (
          storage::Vector< util::Implementation< AssignmentUnaryInterface> >::const_iterator
            itr( GetUnaryAssignments().Begin()), itr_end( GetUnaryAssignments().End());
          itr != itr_end;
          ++itr
        )
        {
          last_instance = util::Enumerated< t_InterfaceType>::AddInstance( new t_ClassType( *itr));
        }
        return last_instance;
      }

    }; // class AssignmentUnaryInterface

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_ASSIGNMENT_UNARY_INTERFACE_H_
