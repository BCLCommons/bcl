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

#ifndef BCL_CLUSTER_NODE_COLORER_H_
#define BCL_CLUSTER_NODE_COLORER_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_cluster_node.h"
#include "util/bcl_util_colors.h"
#include "util/bcl_util_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NodeColorer
    //! @brief NodeColorer takes a node and give the three color points that the node should be
    //! @details The color of the node is determined by calculating a numerical description of the node using some description
    //! function and then converting this numerical description into a color using some coloring function
    //! For example, a node could be colored in a gradient fashion based on the average score of the members the
    //! node contains
    //!
    //! @see @link example_cluster_node_colorer.cpp @endlink
    //! @author alexanns
    //! @date September 6, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class NodeColorer :
      public util::FunctionInterface< Node< t_DataType, t_PrecisionType>, linal::Vector3D>
    {
    protected:

    //////////
    // data //
    //////////

      //! "m_DescriptionFunction" is the function used for getting a numerical description of a Node
      util::ShPtr< util::FunctionInterface< Node< t_DataType, t_PrecisionType>, double> > m_DescriptionFunction;

      //! "m_ColorFunction" is the function used for converting the numerical description of a Node into a color
      //! given by 3 color points
      util::ShPtr< util::FunctionInterface< double, linal::Vector3D> > m_ColorFunction;

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
      NodeColorer() :
        m_DescriptionFunction(),
        m_ColorFunction()
      {
      }

      //! @brief constructor taking parameters for each member variable
      //! @param DESCRIPTION_FUNCTION the function used for getting a numerical description of a Node
      //! @param COLOR_FUNCTION the function used for converting the numerical description of a Node into a color
      //! given by 3 color points
      NodeColorer
      (
        const util::ShPtr< util::FunctionInterface< Node< t_DataType, t_PrecisionType>, double> > &DESCRIPTION_FUNCTION,
        const util::ShPtr< util::FunctionInterface< double, linal::Vector3D> > &COLOR_FUNCTION
      ) :
        m_DescriptionFunction( DESCRIPTION_FUNCTION),
        m_ColorFunction( COLOR_FUNCTION)
      {
      }

      //! @brief virtual copy constructor
      //! @return new copy of this class
      NodeColorer< t_DataType, t_PrecisionType> *Clone() const
      {
        return new NodeColorer< t_DataType, t_PrecisionType>( *this);
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

      //! @brief default color that will be returned if no meaningful numerical description can be determined
      //! @return default color that will be returned if no meaningful numerical description can be determined
      static const util::Color &GetDefaultColor()
      {
        return util::GetColors().e_Black;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking a Triplet and returns a std::string containing the label text
      //! @param NODE the node for which a color will be determined
      //! @return returns the color "NODE" should be
      linal::Vector3D operator()( const Node< t_DataType, t_PrecisionType> &NODE) const
      {
        // true if either "m_ColorFunction" or "m_DescriptionFunction" are not defined ShPtrs
        if( !m_ColorFunction.IsDefined() || !m_DescriptionFunction.IsDefined())
        {
          // return the default color
          return *GetDefaultColor();
        }

        // create const double "numerical_descriptor" and initialize with the value returned by "m_DescriptionFunction"
        const double numerical_descriptor( m_DescriptionFunction->operator()( NODE));

        // true if "numerical_descriptor" is not defined - no descriptors for "NODE" were found
        if( !util::IsDefined( numerical_descriptor))
        {
          // message that default color will be returned
          BCL_MessageStd( "Numerical descriptor undefined.\nReturning default color.");

          // return the default color
          return *GetDefaultColor();
        }

        // return the color points that are determined from "m_DescriptionFunction" and "m_ColorFunction"
        const linal::Vector3D color( m_ColorFunction->operator()( numerical_descriptor));
        BCL_MessageDbg( "color is " + util::Format()( color));
        return color;
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
        io::Serialize::Read( m_DescriptionFunction, ISTREAM);
        io::Serialize::Read( m_ColorFunction, ISTREAM);

        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_DescriptionFunction, OSTREAM, INDENT);
        io::Serialize::Write( m_ColorFunction, OSTREAM, INDENT);

        return OSTREAM;
      }

    }; // template class NodeColorer

    // instantiate s_Instance
    template< typename t_DataType, typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> NodeColorer< t_DataType, t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new NodeColorer< t_DataType, t_PrecisionType>())
    );

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_NODE_COLORER_H_ 
