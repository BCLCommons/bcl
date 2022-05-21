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

#ifndef BCL_PDB_PRINTER_QUALITY_MULTIMER_H_
#define BCL_PDB_PRINTER_QUALITY_MULTIMER_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "quality/bcl_quality_measures.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrinterQualityMultimer
    //! @brief Prints multimer quality measures to PDB lines
    //! @details Takes a protein model, generates the multimer, calculates the given quality measures with respect to
    //!          the multimer native, and returns a table as PDB lines
    //!
    //! @see @link example_pdb_printer_quality_multimer.cpp @endlink
    //! @author weinerbe
    //! @date Nov 17, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PrinterQualityMultimer :
      public util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> >
    {

    private:

    //////////
    // data //
    //////////

      //! set of qualities to be calculated and printed
      storage::Set< quality::Measure> m_Qualities;

      //! native multimer
      util::ShPtr< assemble::ProteinModel> m_NativeMultimer;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PrinterQualityMultimer();

      //! @brief construct from members
      //! @param QUALITIES quality measures to calculate
      //! @param NATIVE native multimer model
      PrinterQualityMultimer
      (
        const storage::Set< quality::Measure> &QUALITIES,
        const util::ShPtr< assemble::ProteinModel> &NATIVE
      );

      //! @brief Clone function
      //! @return pointer to new PrinterQualityMultimer
      PrinterQualityMultimer *Clone() const;

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

    }; // class PrinterQualityMultimer

  } // namespace pdb
} // namespace bcl

#endif // BCL_PDB_PRINTER_QUALITY_MULTIMER_H_ 
