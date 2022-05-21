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
#include "linal/bcl_linal_operations.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_interface.h"
#include "util/bcl_util_enumerate.hpp"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    Operations< t_DataType>::Operations() :
      util::Enumerate< util::ShPtr< OperationsInterface< t_DataType> >, Operations< t_DataType> >( false)
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &Operations< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the default operations class; linal unless -opencl was given with a valid platform
    template< typename t_DataType>
    const OperationsInterface< t_DataType> &Operations< t_DataType>::GetDefaultOperations() const
    {
      return **m_DefaultType;
    }

    //! @brief set the default operations type from given flag
    template< typename t_DataType>
    const typename Operations< t_DataType>::EnumType &Operations< t_DataType>::SetDefaultOperationsType( const EnumType &DEFAULT)
    {
      m_DefaultType = DEFAULT;
      return m_DefaultType;
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Operations< float>;
    template class BCL_API Operations< double>;

  } // namespace linal

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< linal::OperationsInterface<  float> >, linal::Operations<  float> >;
    template class BCL_API Enumerate< ShPtr< linal::OperationsInterface< double> >, linal::Operations< double> >;

  } // namespace util
} // namespace bcl
