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
#include "opencl/bcl_opencl_extensions.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Extensions::Extensions() :
      e_khr_3d_image_writes              ( AddExtension( "cl_khr_3d_image_writes"              )),
      e_khr_byte_addressable_store       ( AddExtension( "cl_khr_byte_addressable_store"       )),
      e_khr_d3d9_sharing                 ( AddExtension( "cl_khr_d3d9_sharing"                 )),
      e_khr_d3d10_sharing                ( AddExtension( "cl_khr_d3d10_sharing"                )),
      e_khr_d3d11_sharing                ( AddExtension( "cl_khr_d3d11_sharing"                )),
      e_khr_icd                          ( AddExtension( "cl_khr_icd"                          )),
      e_khr_fp16                         ( AddExtension( "cl_khr_fp16"                         )),
      e_khr_fp64                         ( AddExtension( "cl_khr_fp64"                         )),
      e_khr_gl_event                     ( AddExtension( "cl_khr_gl_event"                     )),
      e_khr_gl_sharing                   ( AddExtension( "cl_khr_gl_sharing"                   )),
      e_khr_global_int32_base_atomics    ( AddExtension( "cl_khr_global_int32_base_atomics"    )),
      e_khr_global_int32_extended_atomics( AddExtension( "cl_khr_global_int32_extended_atomics")),
      e_khr_int64_base_atomics           ( AddExtension( "cl_khr_int64_base_atomics"           )),
      e_khr_int64_extended_atomics       ( AddExtension( "cl_khr_int64_extended_atomics"       )),
      e_khr_local_int32_base_atomics     ( AddExtension( "cl_khr_local_int32_base_atomics"     )),
      e_khr_local_int32_extended_atomics ( AddExtension( "cl_khr_local_int32_extended_atomics" )),
      e_khr_select_fprounding_mode       ( AddExtension( "cl_khr_select_fprounding_mode"       )),
      e_ext_device_fission               ( AddExtension( "cl_ext_device_fission"               )),
      e_amd_device_attribute_query       ( AddExtension( "cl_amd_device_attribute_query"       )),
      e_amd_event_callback               ( AddExtension( "cl_amd_event_callback"               )),
      e_amd_fp64                         ( AddExtension( "cl_amd_fp64"                         )),
      e_amd_media_ops                    ( AddExtension( "cl_amd_media_ops"                    )),
      e_amd_popcnt                       ( AddExtension( "cl_amd_popcnt"                       )),
      e_amd_printf                       ( AddExtension( "cl_amd_printf"                       )),
      e_amd_vec3                         ( AddExtension( "cl_amd_vec3"                         )),
      e_apple_gl_sharing                 ( AddExtension( "cl_apple_gl_sharing"                 )),
      e_nv_compiler_options              ( AddExtension( "cl_nv_compiler_options"              )),
      e_nv_d3d9_sharing                  ( AddExtension( "cl_nv_d3d9_sharing"                  )),
      e_nv_d3d10_sharing                 ( AddExtension( "cl_nv_d3d10_sharing"                 )),
      e_nv_d3d11_sharing                 ( AddExtension( "cl_nv_d3d11_sharing"                 )),
      e_nv_device_attribute_query        ( AddExtension( "cl_nv_device_attribute_query"        )),
      e_nv_pragma_unroll                 ( AddExtension( "cl_nv_pragma_unroll"                 )),
      e_intel_printf                     ( AddExtension( "cl_intel_printf"                     )),
      e_intel_immediate_execution        ( AddExtension( "cl_intel_immediate_execution"        )),
      e_intel_exec_by_local_thread       ( AddExtension( "cl_intel_exec_by_local_thread"       ))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Extensions::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add extension enum
    //! @param EXTENSION_STRING string for extension, e.g. "cl_amd_fp64"
    //! @return the enum that was added, e_Undefined if there was an error
    Extension Extensions::AddExtension( const std::string &EXTENSION_STRING)
    {
      return AddEnum( EXTENSION_STRING, ExtensionData( EXTENSION_STRING));
    }

    //! @brief get set of extensions from a string, that was queried CL_PLATFORM_EXTENSIONS or CL_DEVICE_EXTENSIONS
    //! @param EXTENSIONS_STRING string containing multiple extensions, whitespace separated
    //! @return set of extensions
    storage::Set< Extension> Extensions::ExtensionsFromString( const std::string &EXTENSIONS_STRING) const
    {
      // all extension strings
      const storage::Vector< std::string> extension_strings( util::SplitString( EXTENSIONS_STRING));
      const storage::Vector< Extension> extension_vector( extension_strings.Begin(), extension_strings.End());
      // set of extensions
      storage::Set< Extension> extensions;
      extensions.InsertElements( extension_vector.Begin(), extension_vector.End());

      // there will be an undefined entry if one of the extension strings is unknown
      if( extensions.Find( e_Undefined) != extensions.End())
      {
        BCL_MessageDbg( "found at least one unknown extension: " + EXTENSIONS_STRING);
        extensions.Erase( e_Undefined);
      }

      // end
      return extensions;
    }

    //! @brief convert set of extensions into one string for output
    //! @param EXTENSIONS set of extensions
    //! @return string like it would result form query CL_PLATFORM_EXTENSIONS
    std::string Extensions::ExtensionsToString( const storage::Set< Extension> &EXTENSIONS)
    {
      std::string descriptor;

      // iterate over all extensions
      for
      (
        storage::Set< Extension>::const_iterator itr( EXTENSIONS.Begin()), itr_end( EXTENSIONS.End());
        itr != itr_end;
        ++itr
      )
      {
        descriptor += itr->GetName() + " ";
      }

      // end
      return descriptor;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief construct on access function for all Extensions
    //! @return reference to only instances of Extensions
    const Extensions &GetExtensions()
    {
      return Extensions::GetEnums();
    }

  } // namespace opencl

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< opencl::ExtensionData, opencl::Extensions>;

  } // namespace util
} // namespace bcl
