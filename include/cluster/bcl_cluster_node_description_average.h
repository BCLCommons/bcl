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

#ifndef BCL_CLUSTER_NODE_DESCRIPTION_AVERAGE_H_
#define BCL_CLUSTER_NODE_DESCRIPTION_AVERAGE_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_cluster_node.h"
#include "util/bcl_util_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NodeDescriptionAverage
    //! @brief NodeDescriptionAverage is for getting a numerical description of a node
    //! @details For this class, the numerical description is the average of all the numerical descriptions of the members
    //! of the node
    //!
    //! @see @link example_cluster_node_description_average.cpp @endlink
    //! @author alexanns
    //! @date September 6, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class NodeDescriptionAverage :
      public util::FunctionInterface< Node< t_DataType, t_PrecisionType>, double>
    {
    protected:

    //////////
    // data //
    //////////

      // "m_DescriptionFunction" is the method which gives the numerical description of one of the members of the node
      util::ShPtr< util::FunctionInterface< t_DataType, double> > m_DescriptionFunction;

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
      NodeDescriptionAverage() :
        m_DescriptionFunction()
      {
      }

      //! @brief constructor taking parameters for each member variable
      //! @param DESCRIPTION_FUNCTION method which gives the numerical description of one of the members of the node
      NodeDescriptionAverage
      (
        const util::ShPtr< util::FunctionInterface< t_DataType, double> > &DESCRIPTION_FUNCTION
      ) :
        m_DescriptionFunction( DESCRIPTION_FUNCTION)
      {
      }

      //! @brief virtual copy constructor
      //! @return new copy of this class
      NodeDescriptionAverage< t_DataType, t_PrecisionType> *Clone() const
      {
        return new NodeDescriptionAverage< t_DataType, t_PrecisionType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking a Node and gives a numerical description of it which is the average of
      //!        all of the numerical descriptions of the members of "NODE"
      //! @param NODE the node for which a numerical description will be calculated
      //! @return returns a double which is the numerical description of "NODE" i.e. the average of all of the
      //! numerical descriptions of the members of "NODE"
      double operator()( const Node< t_DataType, t_PrecisionType> &NODE) const
      {
        // create double "sum" and initialize to zero and it will hold the sum of all the numerical descriptions of
        // the members of "NODE"
        double sum( 0.0);

        // create size_t "number_descriptors" which will count the number of descriptors which goes into "sum"
        size_t number_descriptors( 0);

        // iterate through all the members of "NODE" in order to sum up their numerical descriptions
        for
        (
          typename util::SiPtrList< t_DataType>::const_iterator
            member_itr( NODE.GetMembers().Begin()), member_itr_end( NODE.GetMembers().End());
          member_itr != member_itr_end;
          ++member_itr
        )
        {
          // create const double "" and initialize with the descriptor for the member currently denoted by "member_itr"
          const double current_descriptor( m_DescriptionFunction->operator()( **member_itr));

          // true if "current_descriptor" is defined - should add it to "sum"
          if( util::IsDefined( current_descriptor))
          {
            // add the numerical description of the member currently denoted by "member_itr" to "sum"
            sum += current_descriptor;

            // increase count of "number_descriptors"
            ++number_descriptors;
          }
        }

        // true if no descriptors were able to be found
        if( sum == 0.0 || number_descriptors == 0)
        {
          return util::GetUndefined< double>();
        }

        // return the average of the numerical descriptions of the members of "NODE"
        const double average_description( sum / number_descriptors);
        BCL_MessageStd
        (
          "Average description for node of girth " + util::Format()( NODE.GetGirth()) + " and size " +
          util::Format()( NODE.GetMembers().GetSize()) + " is " + util::Format()( average_description)
        );
        return average_description;
      }

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_DescriptionFunction, ISTREAM);

        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_DescriptionFunction, OSTREAM, INDENT);

        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // NodeDescriptionAverage

    // instantiate s_Instance
    template< typename t_DataType, typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> NodeDescriptionAverage< t_DataType, t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new NodeDescriptionAverage< t_DataType, t_PrecisionType>())
    );

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_NODE_DESCRIPTION_AVERAGE_H_ 
