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

#ifndef BCL_MATH_ASSIGNMENTS_H_
#define BCL_MATH_ASSIGNMENTS_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_assignment_operation_interface.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Assignments
    //! @brief collects all different assignment types, so that they can be iterated through safely
    //!
    //! @see @link example_math_assignments.cpp @endlink
    //! @author mendenjl
    //! @date Feb 06, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Assignments
    {

    public:

      typedef typename storage::Vector
      <
        util::Implementation< AssignmentOperationInterface< t_DataType> >
      >::const_iterator const_iterator;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! undefined default constructor because no one should be constructing this class
      Assignments();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief get all the different assignments
      //! @return vector of all different assignment types
      static const storage::Vector< util::Implementation< AssignmentOperationInterface< t_DataType> > > &GetAssignments();

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief function to add all instances to the enumerated interface for
      //! @tparam t_InterfaceType the enumerated type to add all newly-created t_ClassTypes to
      //! @tparam t_ClassType the actual class instance; must take a single assignment in its constructor
      template< typename t_InterfaceType, typename t_ClassType>
      static util::ObjectInterface *AddInstances()
      {
        // keep a pointer to the last created instance
        util::ObjectInterface *last_instance( NULL);
        for( const_iterator itr( GetAssignments().Begin()), itr_end( GetAssignments().End()); itr != itr_end; ++itr)
        {
          last_instance = util::Enumerated< t_InterfaceType>::AddInstance( new t_ClassType( *itr));
        }
        return last_instance;
      }

    }; // template class Assignments

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Assignments< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Assignments< float>;

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_ASSIGNMENTS_H_
