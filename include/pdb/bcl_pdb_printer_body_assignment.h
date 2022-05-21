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

#ifndef BCL_PDB_PRINTER_BODY_ASSIGNMENT_H_
#define BCL_PDB_PRINTER_BODY_ASSIGNMENT_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_geometry_interface.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrinterBodyAssignment
    //! @brief class prints assignment info to pdb
    //! @details This class determines the assignment information (consisting of which body (density rod) is occupied by which
    //! sse from the pool / protein model and the sse's orientation in the density rod compared to how it should be in
    //! the native structure). 1 equals correct orientation, 0 equals antiparallel orientation. In the output the
    //! density rods are numbered according to the order they are given in the restraints file. In the output the sses
    //! are numbered according to the order in the pool / protein model (i.e. by sequence).
    //!
    //! @see @link example_pdb_printer_body_assignment.cpp @endlink
    //! @author linders, alexanns, karakam, weinerbe
    //! @date Nov 21, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PrinterBodyAssignment :
      public util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> >
    {

    private:

    //////////
    // data //
    //////////

      //! m_Restraints holds the restraints (density map)
      util::ShPtr
      <
        util::ShPtrVector< restraint::Body>
      > m_Restraints;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! remark number used in PDB lines
      static const size_t s_RemarkNumber = 30;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PrinterBodyAssignment();

      //! @brief construct from restraint information
      //! @param RESTRAINTS holds the restraints (density map)
      PrinterBodyAssignment
      (
        const util::ShPtr
        <
          util::ShPtrVector< restraint::Body>
        > &RESTRAINTS
      );

      //! @brief Clone function
      //! @return pointer to new PrinterBodyAssignment
      PrinterBodyAssignment *Clone() const;

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief write relative rotation and translation information between every SSE in the protein model and the body
      //! it is assigned to
      //! @param PROTEIN_MODEL protein model for which the translation and rotation agreement information is determined
      //! @return pdb lines which are written to
      util::ShPtrList< Line> WriteBodySSEAgreementInformation( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief generate assignments from a protein model
      //! @param PROTEIN_MODEL protein model of interest
      //! @return assignments
      storage::Vector< restraint::SSEAssignment>
      GenerateAssignments( const assemble::ProteinModel &PROTEIN_MODEL) const;

    }; // class PrinterBodyAssignment

  } // namespace pdb
} // namespace bcl

#endif // BCL_PDB_PRINTER_BODY_ASSIGNMENT_H_ 
