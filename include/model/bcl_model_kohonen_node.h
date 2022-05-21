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

#ifndef BCL_MODEL_KOHONEN_NODE_H_
#define BCL_MODEL_KOHONEN_NODE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_running_average.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class KohonenNode
    //! @brief is a reference vector in the network. Represents the basic unit the Kohonen network is
    //!        built from.
    //!
    //! @see @link example_model_kohonen_node.cpp @endlink
    //! @author lemmonwa, mendenjl
    //! @date 7/19/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API KohonenNode :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

    private:

      //! the nodes position in the N-dimensions network
      linal::Vector< float> m_Position;

      //! average feature vector
      math::RunningAverage< linal::Vector< float> > m_FeatureAverage;

      //! average result vector
      math::RunningAverage< linal::Vector< float> > m_ResultAverage;

    public:

      // a single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      KohonenNode()
      {
      }

      //! @brief constructor from parameters
      //! @param POSITION the location of this node within the network.
      //! @param FEATURE_VECTOR the initial feature vector for this node
      //! @param RESULT_VECTOR the initial result vector for this node
      KohonenNode
      (
        const linal::Vector< float> &POSITION,
        const linal::Vector< float> &FEATURE_VECTOR,
        const linal::Vector< float> &RESULT_VECTOR
      );

      //! @brief Clone function
      //! @return pointer to new KohonenNode
      KohonenNode *Clone() const
      {
        return new KohonenNode( *this);
      }

    /////////////////
    // data access //
    /////////////////

    public:

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief gets the reference vector for this node
      //! @return the node's reference vector
      const linal::Vector< float> &GetFeatureVector() const
      {
        return m_FeatureAverage.GetAverage();
      }

      //! @brief gets the reference vector for this node
      //! @return the node's reference vector
      const linal::Vector< float> &GetResultVector() const
      {
        return m_ResultAverage.GetAverage();
      }

      //! @brief return the weight associated with the features (= the sum of alphas passed to MapData)
      //! @return the weight associated with the features (= the sum of alphas passed to MapData)
      size_t GetWeight() const
      {
        return m_FeatureAverage.GetWeight();
      }

      //! @brief gets the position of the node within the network
      //! @return the nodes position
      const linal::Vector< float> &GetPosition() const
      {
        return m_Position;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief associates the data item with this node.
      //! @param FEATURE the feature to associate with this node
      //! @param RESULT the result to average with this node
      //! @param ALPHA the current alpha (learning rate)
      void MapData
      (
        const linal::VectorConstInterface< float> &FEATURE,
        const linal::VectorConstInterface< float> &RESULT,
        const double &ALPHA = 1.0
      )
      {
        m_FeatureAverage.AddWeightedObservation( FEATURE, ALPHA);
        m_ResultAverage.AddWeightedObservation( RESULT, ALPHA);
      }

      //! @brief Reset this node
      void Reset()
      {
        m_FeatureAverage.Reset();
        m_ResultAverage.Reset();
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief adds the data from another node to the current node.
      //! @param OTHER the node to add from
      //! @return a reference to this node
      KohonenNode &operator +=( const KohonenNode &OTHER);

      //! @brief adds the data from another node to the current node with a specified weight
      //! @param OTHER the node to add from
      //! @param WEIGHT the weight to given the reference vector from the newly added node
      //! @return a reference to this node
      KohonenNode &AddNodeWithWeight( const KohonenNode &OTHER, const float &WEIGHT);

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

    }; // class KohonenNode

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_KOHONEN_NODE_H_
