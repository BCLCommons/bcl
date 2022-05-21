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

#ifndef BCL_MC_PRINTER_H_
#define BCL_MC_PRINTER_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_printer_default.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Printer
    //! @brief Printer to be used as standard printer in the mc minimization framework
    //!
    //! @see @link example_mc_printer.cpp @endlink
    //! @author fischea, mendenjl
    //! @date Jan 24, 2013
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class Printer :
      public opti::PrinterDefault< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! message level at which to print the step status
      util::Message::MessageLevelEnum m_StepStatusMessageLevel;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Printer() :
        opti::PrinterDefault< t_ArgumentType, t_ResultType>( false, true),
        m_StepStatusMessageLevel( util::Message::e_Standard)
      {
      }

      //! @brief construct from message level and step count
      //! @param STEP_COUNT step count that tracks the mc steps in the approximation process
      //! @param PRINT_ARGUMENT whether to always print the argument
      //! @param PRINT_RESULT whether to always print the result
      //! @param STEP_STATUS_MESSAGE_LEVEL message level at which to print step status
      Printer
      (
        const bool &PRINT_ARGUMENT = false,
        const bool &PRINT_RESULT = true,
        const util::Message::MessageLevel &STEP_STATUS_MESSAGE_LEVEL = util::Message::e_Verbose
      ) :
        opti::PrinterDefault< t_ArgumentType, t_ResultType>( PRINT_ARGUMENT, PRINT_RESULT),
        m_StepStatusMessageLevel( STEP_STATUS_MESSAGE_LEVEL)
      {
      }

      //! @brief Clone function
      //! @return pointer to a new Printer< t_ArgumentType, t_ResultType>
      Printer< t_ArgumentType, t_ResultType> *Clone() const
      {
        return new Printer< t_ArgumentType, t_ResultType>( *this);
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

      //! @brief returns the class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "ResultStepStatus");
        return s_alias;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief prints information concerning the approximation process based on the status of the given tracker
      //! @param TRACKER holds the status of the approximation process
      void Print( const opti::Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        BCL_Message( m_StepStatusMessageLevel, "step status: " + TRACKER.GetStatusLastStep().GetString());

        opti::PrinterDefault< t_ArgumentType, t_ResultType>::Print( TRACKER);
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
        io::Serialize::Read( m_StepStatusMessageLevel, ISTREAM);

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
        io::Serialize::Write( m_StepStatusMessageLevel, OSTREAM, INDENT) << '\n';

        // end
        return OSTREAM;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer( opti::PrinterDefault< t_ArgumentType, t_ResultType>::GetSerializer());
        serializer.SetClassDescription( "Prints any or all of approximations, results, and step status to the screen");
        serializer.AddInitializer
        (
          "step status level",
          "Message level for step status",
          io::Serialization::GetAgent( &m_StepStatusMessageLevel)
        );
        return serializer;
      }

    }; // template class Printer< t_ArgumentType, t_ResultType>

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_PRINTER_H_
