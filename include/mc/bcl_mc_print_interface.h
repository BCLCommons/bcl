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

#ifndef BCL_MC_PRINT_INTERFACE_H_
#define BCL_MC_PRINT_INTERFACE_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_print_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrintInterface
    //! @brief is a class, for writing argument and result for a mc iterative optimization
    //!
    //! @tparam t_ArgumentType argument to the mc optimization
    //! @tparam t_ResultType the type of the objective to optimize on
    //!
    //! @remarks example unnecessary
    //! @author karakam, fischea
    //! @date Feb 23, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class PrintInterface :
      public opti::PrintInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    public:

      //! @brief singleton function to get format object for round number
      //! @return format object for round number
      static const util::Format &GetRoundNumberFormat()
      {
        // initialize static const format eg. "0203" for round #203
        static const util::Format s_format( util::Format().W( 4).Fill( '0').R());

        // end
        return s_format;
      }

      //! @brief singleton function to get format object for stage number
      //! @return format object for stage number
      static const util::Format &GetStageNumberFormat()
      {
        // initialize static const format eg. "00003" for round #3
        static const util::Format s_format( util::Format().W( 5).Fill( '0').R());

        // end
        return s_format;
      }

      //! @brief singleton function to get format object for iteration number
      //! @return format object for iteration number
      static const util::Format &GetIterationNumberFormat()
      {
        // initialize static const format eg. "00203" for iteration #203
        static const util::Format s_format( util::Format().W( 5).Fill( '0').R());

        // end
        return s_format;
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new PrintInterface< t_ArgumentType, t_ResultType>
      virtual PrintInterface< t_ArgumentType, t_ResultType> *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief sets the prefix for the output
      //! @param prefix prefix for the output
      virtual void SetPrefix( const std::string &PREFIX) = 0;

      //! @brief reset and initialize the printer
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @return true if initialization was successful
      virtual void Initialize( const size_t &ROUND_NUMBER) = 0;

      //! @brief reset and initialize the printer
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @param STAGE_NUMBER for multiple optimizations, a different stage number will be passed
      //! @return true if initialization was successful
      virtual void Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER) = 0;

    }; // template class PrintInterface< t_ArgumentType, t_ResultType>

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_PRINT_INTERFACE_H_
