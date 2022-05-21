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

#ifndef BCL_MC_PRINTER_DEFAULT_H_
#define BCL_MC_PRINTER_DEFAULT_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_mc_print_interface.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrinterDefault
    //! @brief default printer to be used Monte Carlo minimization classes
    //! @details This class provides a default implementation for mc::Interface
    //!
    //! @see @link example_mc_printer_default.cpp @endlink
    //! @author karakam, fischea
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

      //! round number
      size_t m_RoundNumber;

      //! stage number
      size_t m_StageNumber;

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
      PrinterDefault() :
        m_RoundNumber( 0),
        m_StageNumber( 0)
      {
      }

      //! @brief Clone function
      //! @return pointer to new PrinterDefault< t_ArgumentType, t_ResultType>
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

      //! @brief return prefix
      //! @return prefix
      const std::string &GetPrefix() const
      {
        // initialize static prefix
        static const std::string s_prefix;

        // end
        return s_prefix;
      }

      //! @brief set prefix to given PREFIX
      //! @param PREFIX new prefix
      void SetPrefix( const std::string &PREFIX)
      {
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief reset and initialize the printer
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @return true if initialization was successful
      void Initialize( const size_t &ROUND_NUMBER)
      {
        m_RoundNumber = ROUND_NUMBER;
      }

      //! @brief reset and initialize the printer
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @param STAGE_NUMBER for multiple optimizations, a different stage number will be passed
      //! @return true if initialization was successful
      void Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER)
      {
        m_RoundNumber = ROUND_NUMBER;
        m_StageNumber = STAGE_NUMBER;
      }

      //! @brief prints information concerning the approximation process based on the status of the given tracker
      //! @param TRACKER holds the status of the approximation process
      void Print( const opti::Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        //only write when verbose level was set
        BCL_MessageVrb
        (
          "Round:\t" + util::Format()( m_RoundNumber) + "\tIteration:\t" + util::Format()( TRACKER.GetIteration())
        );

        // write Argument
        BCL_MessageDbg( "Argument\n" + util::Format()( TRACKER.GetCurrent()->First()));

        // write Result
        BCL_MessageDbg( "Result\n" + util::Format()( TRACKER.GetCurrent()->Second()));
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
        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

    }; // template class PrinterDefault

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> PrinterDefault< t_ArgumentType, t_ResultType>::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterDefault< t_ArgumentType, t_ResultType>())
    );

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_PRINTER_DEFAULT_H_
