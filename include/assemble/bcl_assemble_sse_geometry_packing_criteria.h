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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_CRITERIA_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_CRITERIA_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "contact/bcl_contact_types.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackingCriteriaCombine
    //! @brief class allows using multiple criteria to evaluate a SSEGeometryPacking
    //! @details This class allows to evaluate the SSEGeometryPacking by looking at combination of multiple criteria.
    //! All the criteria need to return true for this class to return true, otherwise it returns false
    //!
    //! @see @link example_assemble_sse_geometry_packing_criteria.cpp @endlink
    //! @author karakam
    //! @date Jul 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPackingCriteriaCombine :
      public math::FunctionInterfaceSerializable< SSEGeometryPacking, bool>
    {

    private:

    //////////
    // data //
    //////////

      //! vector of criteria
      util::ShPtrVector< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> > m_CriteriaVector;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEGeometryPackingCriteriaCombine();

      //! @brief constructor from a criteria vector
      //! @param CRITERIA_VECTOR vector of criteria
      SSEGeometryPackingCriteriaCombine
      (
        const util::ShPtrVector< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> > &CRITERIA_VECTOR
      );

      //! @brief virtual copy constructor
      SSEGeometryPackingCriteriaCombine *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate the given SSEGeometryPacking according to the group of criteria
      //! @param SSE_GEOMETRY_PACKING SSEGeometryPacking to be evaluated
      //! @return true if given SSEGeometryPacking follows all criteria
      bool operator()( const SSEGeometryPacking &SSE_GEOMETRY_PACKING) const;

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

    }; // class SSEGeometryPackingCriteriaCombine

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackingCriteriaDistance
    //! @brief class checks if the distance of SSEGeometryPacking is valid
    //! @details This class allows to evaluate the SSEGeometryPacking by looking at whether its' distance meets a cutoff
    //! It can be used with any comparison operator
    //!
    //! @see @link example_assemble_sse_geometry_packing_criteria.cpp @endlink
    //! @author karakam
    //! @date Jul 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPackingCriteriaDistance :
      public math::FunctionInterfaceSerializable< SSEGeometryPacking, bool>
    {

    private:

    //////////
    // data //
    //////////

      //! distance cutoff
      double m_DistanceCutoff;

      //! operator to be used
      math::Comparisons< double>::Comparison m_Comparison;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEGeometryPackingCriteriaDistance();

      //! @brief constructor from a distance cutoff and a comparison
      //! @param DISTANCE_CUTOFF distance cutoff
      //! @param COMPARISON Comparison operator
      SSEGeometryPackingCriteriaDistance
      (
        const double DISTANCE_CUTOFF,
        const math::Comparisons< double>::Comparison &COMPARISON
      );

      //! @brief virtual copy constructor
      SSEGeometryPackingCriteriaDistance *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate the given SSEGeometryPacking according to distance cutoff criteria
      //! @param SSE_GEOMETRY_PACKING SSEGeometryPacking to be evaluated
      //! @return true if SSEGeometryPacking follows distance cutoff criteria
      bool operator()( const SSEGeometryPacking &SSE_GEOMETRY_PACKING) const;

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

    }; // class SSEGeometryPackingCriteriaDistance

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackingCriteriaContactType
    //! @brief class checks if the given SSEGeometryPacking has a valid contact type
    //! @details This class allows to evaluate the SSEGeometryPacking by looking at whether its' contact type is one of the
    //! the specified contact types.
    //!
    //! @see @link example_assemble_sse_geometry_packing_criteria.cpp @endlink
    //! @author karakam
    //! @date Jul 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPackingCriteriaContactType :
      public math::FunctionInterfaceSerializable< SSEGeometryPacking, bool>
    {

    private:

    //////////
    // data //
    //////////

      //! contact types
      storage::Set< contact::Type> m_ContactTypes;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEGeometryPackingCriteriaContactType();

      //! @brief constructor from a contact type
      //! @param CONTACT_TYPE requested contact type
      SSEGeometryPackingCriteriaContactType( const contact::Type &CONTACT_TYPE);

      //! @brief constructor from a set of contact types
      //! @param CONTACT_TYPES set of requested contact types
      SSEGeometryPackingCriteriaContactType( const storage::Set< contact::Type> &CONTACT_TYPES);

      //! @brief virtual copy constructor
      SSEGeometryPackingCriteriaContactType *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate the given SSEGeometryPacking according to contact type criteria
      //! @param SSE_GEOMETRY_PACKING SSEGeometryPacking to be evaluated
      //! @return true if SSEGeometryPacking follows contact type criteria
      bool operator()( const SSEGeometryPacking &SSE_GEOMETRY_PACKING) const;

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

    }; // class SSEGeometryPackingCriteriaContactType

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackingCriteriaInteractionWeight
    //! @brief class checks if the interaction weight of SSEGeometryPacking is valid
    //! @details This class allows to evaluate the SSEGeometryPacking by looking at whether its' interaction weight meets a cutoff
    //! It can be used with any comparison operator
    //!
    //! @see @link example_assemble_sse_geometry_packing_criteria.cpp @endlink
    //! @author karakam
    //! @date Jul 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPackingCriteriaInteractionWeight :
      public math::FunctionInterfaceSerializable< SSEGeometryPacking, bool>
    {

    private:

    //////////
    // data //
    //////////

      //! interaction weight cutoff
      double m_InteractionWeightCutoff;

      //! operator to be used
      math::Comparisons< double>::Comparison m_Comparison;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEGeometryPackingCriteriaInteractionWeight();

      //! @brief constructor from an interaction weight cutoff and a comparison
      //! @param INTERACTION_WEIGHT_CUTOFF interaction weight cutoff
      //! @param COMPARISON Comparison operator
      SSEGeometryPackingCriteriaInteractionWeight
      (
        const double INTERACTION_WEIGHT_CUTOFF,
        const math::Comparisons< double>::Comparison &COMPARISON
      );

      //! @brief virtual copy constructor
      SSEGeometryPackingCriteriaInteractionWeight *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate the given SSEGeometryPacking according to interaction weight cutoff criteria
      //! @param SSE_GEOMETRY_PACKING SSEGeometryPacking to be evaluated
      //! @return true if SSEGeometryPacking follows interaction weight cutoff criteria
      bool operator()( const SSEGeometryPacking &SSE_GEOMETRY_PACKING) const;

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

    }; // class SSEGeometryPackingCriteriaInteractionWeight

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackingCriteriaStrandWeight
    //! @brief class checks if the strand pairing weight of SSEGeometryPacking is valid
    //! @details This class allows to evaluate the SSEGeometryPacking by looking at whether its' strand pairing weight meets a
    //! cutoff.
    //! It can be used with any comparison operator
    //!
    //! @see @link example_assemble_sse_geometry_packing_criteria.cpp @endlink
    //! @author karakam
    //! @date Jul 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPackingCriteriaStrandWeight :
      public math::FunctionInterfaceSerializable< SSEGeometryPacking, bool>
    {

    private:

    //////////
    // data //
    //////////

      //! strand weight cutoff
      double m_StrandWeightCutoff;

      //! operator to be used
      math::Comparisons< double>::Comparison m_Comparison;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEGeometryPackingCriteriaStrandWeight();

      //! @brief constructor from a strand weight cutoff and a comparison
      //! @param STRAND_WEIGHT_CUTOFF strand weight cutoff
      //! @param COMPARISON Comparison operator
      SSEGeometryPackingCriteriaStrandWeight
      (
        const double STRAND_WEIGHT_CUTOFF,
        const math::Comparisons< double>::Comparison &COMPARISON
      );

      //! @brief virtual copy constructor
      SSEGeometryPackingCriteriaStrandWeight *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate the given SSEGeometryPacking according to strand weight cutoff criteria
      //! @param SSE_GEOMETRY_PACKING SSEGeometryPacking to be evaluated
      //! @return true if SSEGeometryPacking follows strand weight cutoff criteria
      bool operator()( const SSEGeometryPacking &SSE_GEOMETRY_PACKING) const;

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

    }; // class SSEGeometryPackingCriteriaStrandWeight

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackingCriteriaDistancePerType
    //! @brief checks to see if the distance is valid considering the contact type
    //! @details compares the distance to see if it is within the acceptable distance range for the contact type
    //!
    //! @see @link example_assemble_sse_geometry_packing_criteria.cpp @endlink
    //! @author weinerbe
    //! @date Sep 17, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPackingCriteriaDistancePerType :
      public math::FunctionInterfaceSerializable< SSEGeometryPacking, bool>
    {

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
      SSEGeometryPackingCriteriaDistancePerType();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief virtual copy constructor
      SSEGeometryPackingCriteriaDistancePerType *Clone() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate the given SSEGeometryPacking according to the distance criteria
      //! @param SSE_GEOMETRY_PACKING SSEGeometryPacking to be evaluated
      //! @return true if SSEGeometryPacking follows the distance criteria
      bool operator()( const SSEGeometryPacking &SSE_GEOMETRY_PACKING) const;

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

    }; // class SSEGeometryPackingCriteriaDistancePerType

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_CRITERIA_H_ 
