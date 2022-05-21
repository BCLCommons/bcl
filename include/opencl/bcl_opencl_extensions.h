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

#ifndef BCL_OPENCL_EXTENSIONS_H_
#define BCL_OPENCL_EXTENSIONS_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_extension_data.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Extensions
    //! @brief an enumerator on extension data objects
    //! @details extenions enumerates all possible opencl extensions for easy access and usability. It also provides
    //!          handling of extension strings, that are white space separted as they can be queried from the platform
    //!          or device
    //!
    //! @see @link example_opencl_extensions.cpp @endlink
    //! @author woetzen
    //! @date Aug 24, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Extensions :
      public util::Enumerate< ExtensionData, Extensions>
    {

      friend class util::Enumerate< ExtensionData, Extensions>;

    public:

    //////////
    // data //
    //////////

      // khronos
      Extension e_khr_3d_image_writes              ;
      Extension e_khr_byte_addressable_store       ;
      Extension e_khr_d3d9_sharing                 ;
      Extension e_khr_d3d10_sharing                ;
      Extension e_khr_d3d11_sharing                ;
      Extension e_khr_icd                          ;
      Extension e_khr_fp16                         ;
      Extension e_khr_fp64                         ;
      Extension e_khr_gl_event                     ;
      Extension e_khr_gl_sharing                   ;
      Extension e_khr_global_int32_base_atomics    ;
      Extension e_khr_global_int32_extended_atomics;
      Extension e_khr_int64_base_atomics           ;
      Extension e_khr_int64_extended_atomics       ;
      Extension e_khr_local_int32_base_atomics     ;
      Extension e_khr_local_int32_extended_atomics;
      Extension e_khr_select_fprounding_mode       ;

      // other extensions
      Extension e_ext_device_fission               ;

      // amd specific
      Extension e_amd_device_attribute_query       ;
      Extension e_amd_event_callback               ;
      Extension e_amd_fp64                         ;
      Extension e_amd_media_ops                    ;
      Extension e_amd_popcnt                       ;
      Extension e_amd_printf                       ;
      Extension e_amd_vec3                         ;

      // apple specific
      Extension e_apple_gl_sharing                 ;

      // nvidia specific
      Extension e_nv_compiler_options              ;
      Extension e_nv_d3d9_sharing                  ;
      Extension e_nv_d3d10_sharing                 ;
      Extension e_nv_d3d11_sharing                 ;
      Extension e_nv_device_attribute_query        ;
      Extension e_nv_pragma_unroll                 ;

      // intel specific
      Extension e_intel_printf                     ;
      Extension e_intel_immediate_execution        ;
      Extension e_intel_exec_by_local_thread       ;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Extensions();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief add extension enum
      //! @param EXTENSION_STRING string for extension, e.g. "cl_amd_fp64"
      //! @return the enum that was added, e_Undefined if there was an error
      Extension AddExtension( const std::string &EXTENSION_STRING);

      //! @brief get set of extensions from a string, that was queried CL_PLATFORM_EXTENSIONS or CL_DEVICE_EXTENSIONS
      //! @param EXTENSIONS_STRING string containing multiple extensions, whitespace separated
      //! @return set of extensions
      storage::Set< Extension> ExtensionsFromString( const std::string &EXTENSIONS_STRING) const;

      //! @brief convert set of extensions into one string for output
      //! @param EXTENSIONS set of extensions
      //! @return string like it would result form query CL_PLATFORM_EXTENSIONS
      static std::string ExtensionsToString( const storage::Set< Extension> &EXTENSIONS);

    }; // class Extensions

    //! @brief construct on access function for all Extensions
    //! @return reference to only instances of Extensions
    BCL_API
    const Extensions &GetExtensions();

  } // namespace opencl

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< opencl::ExtensionData, opencl::Extensions>;

  } // namespace util
} // namespace bcl

#endif // BCL_OPENCL_EXTENSIONS_H_ 
