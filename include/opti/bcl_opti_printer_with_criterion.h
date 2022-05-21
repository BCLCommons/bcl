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

#ifndef BCL_OPTI_PRINTER_WITH_CRITERION_H_
#define BCL_OPTI_PRINTER_WITH_CRITERION_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_criterion_interface.h"
#include "bcl_opti_print_interface.h"
#include "bcl_opti_printer_default.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrinterWithCriterion
    //! @brief Conditional printing
    //!
    //! @see @link example_opti_printer_with_criterion.cpp @endlink
    //! @author mendenjl
    //! @date Sep 03, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class PrinterWithCriterion :
      public PrintInterface< t_ArgumentType, t_ResultType>
    {

    private:

    //////////
    // data //
    //////////

      //! base printer to be called
      util::Implementation< PrintInterface< t_ArgumentType, t_ResultType> > m_Printer;

      //! test whether to call the printer
      util::Implementation< CriterionInterface< t_ArgumentType, t_ResultType> > m_PrintCriteria;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PrinterWithCriterion() :
        m_Printer( new PrinterDefault< t_ArgumentType, t_ResultType>())
      {
      }

      //! @brief constructor from message levels for argument and result
      //! @param PRINTER true to print the argument
      //! @param PRINT_RESULT true to print the result
      PrinterWithCriterion
      (
        const util::Implementation< PrintInterface< t_ArgumentType, t_ResultType> > &PRINTER,
        const util::Implementation< CriterionInterface< t_ArgumentType, t_ResultType> > &PRINT_CRITERIA
      ) :
        m_Printer( PRINTER),
        m_PrintCriteria( PRINT_CRITERIA)
      {
      }

      //! @brief Clone function
      //! @return pointer to new PrinterWithCriterion< t_ArgumentType, typename t_ResultType>
      PrinterWithCriterion< t_ArgumentType, t_ResultType> *Clone() const
      {
        return new PrinterWithCriterion< t_ArgumentType, t_ResultType>( *this);
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
        static const std::string s_alias( "Conditional");
        return s_alias;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief prints information concerning the approximation process based on the status of the given tracker
      //! @param TRACKER holds the status of the approximation process
      void Print( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        if( !m_PrintCriteria.IsDefined() || m_PrintCriteria->CriteriaMet( TRACKER))
        {
          m_Printer->Print( TRACKER);
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
        serializer.SetClassDescription( "Allows a printer to only be called when selected criteria are met");
        serializer.AddInitializer
        (
          "printer",
          "printer to call",
          io::Serialization::GetAgent( &m_Printer)
        );
        serializer.AddInitializer
        (
          "condition",
          "Criteria that tells when to call the printer",
          io::Serialization::GetAgent( &m_PrintCriteria)
        );
        return serializer;
      }

    }; // template class PrinterWithCriterion

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> PrinterWithCriterion< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< PrintInterface< t_ArgumentType, t_ResultType> >::AddInstance( new PrinterWithCriterion< t_ArgumentType, t_ResultType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_PRINTER_WITH_CRITERION_H_
