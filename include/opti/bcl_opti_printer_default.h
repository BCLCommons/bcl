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

#ifndef BCL_OPTI_PRINTER_DEFAULT_H_
#define BCL_OPTI_PRINTER_DEFAULT_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_print_interface.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrinterDefault
    //! @brief default printer to be used in optimization framework
    //! @details default PrinterInterface implementation that prints the argument and results to the output stream
    //!
    //! @see @link example_opti_printer_default.cpp @endlink
    //! @author karakam
    //! @date Feb 23, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class PrinterDefault :
      public PrintInterface< t_ArgumentType, t_ResultType>
    {

    private:

    //////////
    // data //
    //////////

      bool m_PrintArgument;  //!< Whether to print the argument
      bool m_PrintResult;    //!< Whether to print the result

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PrinterDefault() :
        m_PrintArgument( false),
        m_PrintResult( false)
      {
      }

      //! @brief constructor from message levels for argument and result
      //! @param PRINT_ARGUMENT true to print the argument
      //! @param PRINT_RESULT true to print the result
      PrinterDefault
      (
        const bool &PRINT_ARGUMENT,
        const bool &PRINT_RESULT
      ) :
        m_PrintArgument( PRINT_ARGUMENT),
        m_PrintResult( PRINT_RESULT)
      {
      }

      //! @brief Clone function
      //! @return pointer to new PrinterDefault< t_ArgumentType, typename t_ResultType>
      PrinterDefault< t_ArgumentType, t_ResultType> *Clone() const
      {
        return new PrinterDefault< t_ArgumentType, t_ResultType>( *this);
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
        static const std::string s_alias( "Display");
        return s_alias;
      }

      //! @brief get whether the result is to be display
      //! @return true if the result should be displayed
      const bool &GetDisplayResult() const
      {
        return m_PrintResult;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief prints information concerning the approximation process based on the status of the given tracker
      //! @param TRACKER holds the status of the approximation process
      void Print( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        if( TRACKER.GetPhase() != e_Iteration)
        {
          // tracker default has no use for the start model
          return;
        }
        if( m_PrintResult)
        {
          // write iteration
          BCL_MessageCrt
          (
            "Iteration: " + util::Format()( TRACKER.GetIteration())
            + " " + TRACKER.GetStatusOfLastStep().GetString()
            + " Result: " + util::Format()( TRACKER.GetCurrent()->Second())
            + " Best: "   + util::Format()( TRACKER.GetBest()->Second())
          );
        }
        if( m_PrintArgument)
        {
          // write Argument
          BCL_MessageCrt( "Argument\n" + util::Format()( TRACKER.GetCurrent()->First()));
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Prints arguments and / or results to the screen");
        serializer.AddInitializer
        (
          "approximation",
          "True to print the approximation to the screen at each iteration",
          io::Serialization::GetAgent( &m_PrintArgument),
          "False"
        );
        serializer.AddInitializer
        (
          "result",
          "True to print the score/result/objective function at each iteration",
          io::Serialization::GetAgent( &m_PrintResult),
          "True"
        );
        return serializer;
      }

    }; // template class PrinterDefault

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> PrinterDefault< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< PrintInterface< t_ArgumentType, t_ResultType> >::AddInstance( new PrinterDefault< t_ArgumentType, t_ResultType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_PRINTER_DEFAULT_H_
