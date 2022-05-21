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

#ifndef BCL_ALIGN_HANDLER_CLASSES_H_
#define BCL_ALIGN_HANDLER_CLASSES_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_handler_blc.h"
#include "bcl_align_handler_fasta.h"
#include "bcl_align_handler_interface.h"
#include "bcl_align_handler_pir.h"
#include "bcl_align_handler_standard.h"
#include "biol/bcl_biol_aa_base.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerClasses
    //! @brief enumerates all alignment handler
    //!
    //! @see @link example_align_handler_classes.cpp @endlink
    //! @author heinzes1
    //! @date Jan 25, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class HandlerClasses :
      public util::Enumerate< util::ShPtr< HandlerInterface< t_Member> >, HandlerClasses< t_Member> >
    {
      friend class util::Enumerate< util::ShPtr< HandlerInterface< t_Member> >, HandlerClasses< t_Member> >;

    public:

    //////////
    // data //
    //////////

      //! @brief for HandlerClasses
      typedef util::Enum< util::ShPtr< HandlerInterface< t_Member> >, HandlerClasses< t_Member> > HandlerClass;

      // declare all handler classes
      const HandlerClass e_BLC;
      const HandlerClass e_Fasta;
      const HandlerClass e_PIR;
      const HandlerClass e_Standard;

      //! @brief gives the flag which allows specifying output formats that are desired
      //! @return flag which allows specifying output formats that are desired
      static const util::ShPtr< command::FlagInterface> &GetFlagOutputFormats();

      //! @brief gives the flag which allows specifying the input format for an alignment
      //! @return flag which allows specifying the input format for an alignment
      static const util::ShPtr< command::FlagInterface> &GetFlagInputFormat();

      using util::Enumerate< util::ShPtr< HandlerInterface< t_Member> >, HandlerClasses< t_Member> >::AddEnum;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all HandlerClasses
      HandlerClasses() :
        e_BLC(      AddEnum( "blc"     , util::ShPtr< HandlerInterface< t_Member> >( new HandlerBLC< t_Member>()))),
        e_Fasta(    AddEnum( "fasta"   , util::ShPtr< HandlerInterface< t_Member> >( new HandlerFasta< t_Member>()))),
        e_PIR(      AddEnum( "pir"     , util::ShPtr< HandlerInterface< t_Member> >( new HandlerPIR< t_Member>()))),
        e_Standard( AddEnum( "standard", util::ShPtr< HandlerInterface< t_Member> >( new HandlerStandard< t_Member>( util::ShPtr< score::AssignmentWithGap< t_Member> >()))))
      {
      }

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    }; // template class HandlerClasses

    //! @brief construct access function for all HandlerClasses
    //! @return reference to instances of HandlerClasses
    template< typename t_Member>
    const HandlerClasses< t_Member> &GetHandlerClasses()
    {
      return HandlerClasses< t_Member>::GetEnums();
    }

    //! @brief gives the flag which allows specifying output formats that are desired
    //! @return flag which allows specifying output formats that are desired
    template< typename t_Member>
    const util::ShPtr< command::FlagInterface> &HandlerClasses< t_Member>::GetFlagOutputFormats()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "outputformat",
          "output formats the alignment should be printed in",
          command::Parameter
          (
            "outputformat_type",
            "output formats",
            command::ParameterCheckEnumerate< HandlerClasses< biol::AABase> >(),
            GetHandlerClasses< biol::AABase>().e_PIR.GetName()
          ),
          0,
          GetHandlerClasses< biol::AABase>().GetEnumCount()
        )
      );

      return s_flag;
    }

    //! @brief gives the flag which allows specifying the input format for an alignment
    //! @return flag which allows specifying the input format for an alignment
    template< typename t_Member>
    const util::ShPtr< command::FlagInterface> &HandlerClasses< t_Member>::GetFlagInputFormat()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "inputformat",
          "input format to read in the alignment",
          command::Parameter
          (
            "inputformat_type",
            "input format",
            command::ParameterCheckEnumerate< HandlerClasses< biol::AABase> >(),
            GetHandlerClasses< biol::AABase>().e_PIR.GetName()
          )
        )
      );

      return s_flag;
    }

  } // namespace align

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< align::HandlerInterface< biol::AABase> >, align::HandlerClasses< biol::AABase> >;

  } // namespace util
} // namespace bcl

#endif // BCL_ALIGN_HANDLER_CLASSES_H_ 
