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

#ifndef BCL_MC_PRINTER_COMBINED_H_
#define BCL_MC_PRINTER_COMBINED_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_mc_print_interface.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrinterCombined
    //! @brief class that allows to combine several printers
    //! @details This class combines several printers ( e.g. the score printer and the body assignment printer) allowing to print
    //! different things at once.
    //!
    //! @see @link example_mc_printer_combined.cpp @endlink
    //! @author linders, alexanns
    //! @date Feb 10, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class PrinterCombined :
      public PrintInterface< t_ArgumentType, t_ResultType>
    {

    private:

    //////////
    // data //
    //////////

      //! prefix for the combined printer
      std::string m_Prefix;

      //! list of printers that PrinterCombined contains
      util::ShPtrList< PrintInterface< t_ArgumentType, t_ResultType> > m_Printers;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PrinterCombined() :
        m_Printers()
      {
      }

      //! @brief Clone function
      //! @return pointer to new PrinterCombined
      PrinterCombined< t_ArgumentType, t_ResultType> *Clone() const
      {
        return new PrinterCombined< t_ArgumentType, t_ResultType>( *this);
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

      //! @brief access to all printers
      //! @return list of printer interfaces, that are combined
      const util::ShPtrList< PrintInterface< t_ArgumentType, t_ResultType> > &GetPrinters() const
      {
        return m_Printers;
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief return prefix
      //! @return prefix
      const std::string &GetPrefix() const
      {
        return m_Prefix;
      }

      //! @brief set prefix to given PREFIX
      //! @param PREFIX new prefix
      void SetPrefix( const std::string &PREFIX)
      {
        // set prefix for this class
        m_Prefix = PREFIX;

        // iterate over all the printer
        for
        (
          typename util::ShPtrList< PrintInterface< t_ArgumentType, t_ResultType> >::iterator
            printer_itr( m_Printers.Begin()), printer_itr_end( m_Printers.End());
          printer_itr != printer_itr_end; ++printer_itr
        )
        {
          // and set the prefix on this printer
          ( *printer_itr)->SetPrefix( PREFIX);
        }
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief reset and initialize the printer
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @return true if initialization was successful
      void Initialize( const size_t &ROUND_NUMBER)
      {
        // iterate over all the printers that PrinterCombined contains
        for
        (
          typename util::ShPtrList< PrintInterface< t_ArgumentType, t_ResultType> >::iterator
            printer_itr( m_Printers.Begin()), printer_itr_end( m_Printers.End());
          printer_itr != printer_itr_end; ++printer_itr
        )
        {
          // call the initialize function for the individual printers
          ( *printer_itr)->Initialize( ROUND_NUMBER);
        }
      }

      //! @brief reset and initialize the printer
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @param STAGE_NUMBER for multiple stage optimizations, a different round and stage number will be passed
      //! @return true if initialization was successful
      void Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER)
      {
        // iterate over all the printers that PrinterCombined contains
        for
        (
          typename util::ShPtrList< PrintInterface< t_ArgumentType, t_ResultType> >::iterator
            printer_itr( m_Printers.Begin()), printer_itr_end( m_Printers.End());
          printer_itr != printer_itr_end; ++printer_itr
        )
        {
          // call the initialize function for the individual printers
          ( *printer_itr)->Initialize( ROUND_NUMBER, STAGE_NUMBER);
        }
      }

      //! @brief prints information concerning the approximation process based on the status of the given tracker
      //! @param TRACKER holds the status of the approximation process
      void Print( const opti::Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        // iterate over all the printers that PrinterCombined contains
        for
        (
          typename util::ShPtrList< PrintInterface< t_ArgumentType, t_ResultType> >::const_iterator
            printer_itr( m_Printers.Begin()), printer_itr_end( m_Printers.End());
          printer_itr != printer_itr_end; ++printer_itr
        )
        {
          // call the print function for the individual printers
          ( *printer_itr)->Print( TRACKER);
        }
      }

      //! @brief insert a printer into the list of printers that PrinterCombined contains
      //! @param PRINTER ShPtr to printer that is to be inserted
      void Insert( const util::ShPtr< PrintInterface< t_ArgumentType, t_ResultType> > &PRINTER)
      {
        // insert PRINTER into m_Printers
        m_Printers.InsertElement( PRINTER);
      }

      //! @brief remove a printer from the list of printers that PrinterCombined contains
      //! @param PRINTER ShPtr to printer that is to be removed
      void Remove( const util::ShPtr< PrintInterface< t_ArgumentType, t_ResultType> > &PRINTER)
      {
        // get iterator to printer that you want to remove from list using the std::find function
        typename util::ShPtrList< PrintInterface< t_ArgumentType, t_ResultType> >::iterator
          itr_to_element( std::find( m_Printers.Begin(), m_Printers.End(), PRINTER));

        // if found
        if( itr_to_element != m_Printers.End())
        {
          // remove located printer
          m_Printers.RemoveElement( itr_to_element);
        }
        else
        {
          BCL_MessageStd( "Can't find the provided printer to remove")
        }
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
        // read members
        io::Serialize::Read( m_Prefix, ISTREAM);
        io::Serialize::Read( m_Printers, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_Prefix, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Printers, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class PrinterCombined

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> PrinterCombined< t_ArgumentType, t_ResultType>::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterCombined< t_ArgumentType, t_ResultType>())
    );

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_PRINTER_COMBINED_H_
