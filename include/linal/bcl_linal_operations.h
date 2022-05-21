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

#ifndef BCL_LINAL_OPERATIONS_H_
#define BCL_LINAL_OPERATIONS_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_operations_interface.h"
#include "command/bcl_command_flag_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

using bcl::util::ShPtr;

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Operations
    //! @brief enumerates implementations of OperationsInterface
    //!
    //! @see @link example_linal_operations.cpp @endlink
    //! @author woetzen, vuot2
    //! @date Mar 14, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Operations :
      public util::Enumerate< util::ShPtr< OperationsInterface< t_DataType> >, Operations< t_DataType> >
    {
      friend class util::Enumerate< util::ShPtr< OperationsInterface< t_DataType> >, Operations< t_DataType> >;

    public:

      typedef typename util::Enumerate< util::ShPtr< OperationsInterface< t_DataType> >, Operations< t_DataType> >::EnumType Operator;

    private:
      Operator m_DefaultType; //!< Default value; set over command line

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Operations();

    public:

      typedef typename util::Enumerate< util::ShPtr< OperationsInterface< t_DataType> >, Operations< t_DataType> >::EnumType EnumType;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the default operations class; linal unless -opencl was given with a valid platform
      const OperationsInterface< t_DataType> &GetDefaultOperations() const;

      //! @brief set the default operations type from given flag
      //! @return the passed in enum
      const EnumType &SetDefaultOperationsType( const EnumType &DEFAULT);

    }; // class Operations

    template< typename t_DataType>
    inline
    const Operations< t_DataType> &GetOperations()
    {
      return Operations< t_DataType>::GetEnums();
    }

    template< typename t_DataType>
    inline
    Operations< t_DataType> &GetOperationsNonConst()
    {
      return Operations< t_DataType>::GetEnums();
    }

    template< typename t_DataType>
    inline
    const OperationsInterface< t_DataType> &GetDefaultOperations()
    {
      return Operations< t_DataType>::GetEnums().GetDefaultOperations();
    }

    BCL_EXPIMP_TEMPLATE template class BCL_API Operations< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Operations< double>;

  } // namespace linal

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< linal::OperationsInterface<  float> >, linal::Operations<  float> >;
    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< linal::OperationsInterface< double> >, linal::Operations< double> >;

  } // namespace util
} // namespace bcl

#endif // BCL_LINAL_OPERATIONS_H_
