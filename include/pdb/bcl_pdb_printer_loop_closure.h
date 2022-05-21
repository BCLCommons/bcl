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

#ifndef BCL_PDB_PRINTER_LOOP_CLOSURE_H_
#define BCL_PDB_PRINTER_LOOP_CLOSURE_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrinterLoopClosure
    //! @brief Prints information to pdb lines about loop closing
    //! @details Prints to pdb lines the loop closure cutoff, the distance, and whether or not the loop was closed
    //!           according to these two values.
    //!
    //! @see @link example_pdb_printer_loop_closure.cpp @endlink
    //! @author weinerbe
    //! @date Dec 12, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PrinterLoopClosure :
      public util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> >
    {

    private:

    //////////
    // data //
    //////////

      //! the threshold for considering the loop closed
      double m_ClosureThreshold;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! remark number used in PDB lines
      static const size_t s_RemarkNumber = 20;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PrinterLoopClosure();

      //! @brief constructor taking member variable parameters
      //! @param CLOSURE_THRESHOLD the threshold for considering the loop closed
      PrinterLoopClosure( const double CLOSURE_THRESHOLD);

      //! @brief Clone function
      //! @return pointer to new PrinterLoopClosure
      PrinterLoopClosure *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief takes a protein model and returns PDB lines
      //! @param PROTEIN_MODEL protein model to print
      //! @return PDB lines
      util::ShPtrList< Line> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class PrinterLoopClosure

  } // namespace pdb
} // namespace bcl

#endif // BCL_PDB_PRINTER_LOOP_CLOSURE_H_ 
