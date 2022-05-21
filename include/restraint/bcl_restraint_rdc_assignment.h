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

#ifndef BCL_RESTRAINT_RDC_ASSIGNMENT_H_
#define BCL_RESTRAINT_RDC_ASSIGNMENT_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RDCAssignment
    //! @brief Stores a collection of coordinates and the associated experimental value.
    //! @details Contains a collection of coordinates and experimental RDC values.  The coordinates are typically
    //!          obtained from the GenerateAssignment function in the RDC class.  Other classes, such as
    //!          nmr::ResidualDipolarCouplingLeastSquareDeviation can use this data to generate the
    //!          theoretical RDC values.
    //!
    //! @see @link example_restraint_rdc_assignment.cpp @endlink
    //! @author weinerbe
    //! @date Mar 16, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RDCAssignment :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! list of atom coordinates and experimental RDC value
      storage::List< storage::Triplet< linal::Vector3D, linal::Vector3D, double> > m_Data;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RDCAssignment();

      //! @brief Clone function
      //! @return pointer to new RDCAssignment
      RDCAssignment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the RDC assignment data
      //! @return the RDC assignment data
      const storage::List< storage::Triplet< linal::Vector3D, linal::Vector3D, double> > &GetData() const
      {
        return m_Data;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief add new data into the list
      //! @param FIRST_COORD first coordinate
      //! @param SECOND_COORD second coordinate
      //! @param VALUE experimental RDC value
      void PushBack
      (
        const linal::Vector3D &FIRST_COORD,
        const linal::Vector3D &SECOND_COORD,
        const double VALUE
      );

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

    }; // class RDCAssignment

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_RDC_ASSIGNMENT_H_ 
