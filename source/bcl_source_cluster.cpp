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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "cluster/bcl_cluster.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace cluster
} // namespace bcl
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

//// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "cluster/bcl_cluster_input_classes.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {

  } // namespace cluster

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< cluster::InputInterface< std::string, float> >, cluster::InputClasses< std::string, float> >;
    template class BCL_API Enumerate< ShPtr< cluster::InputInterface< linal::Vector< float>, float> >, cluster::InputClasses< linal::Vector< float>, float> >;

  } // namespace util
} // namespace bcl
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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "cluster/bcl_cluster_input_features.h"

// includes from bcl - sorted alphabetically
#include "cluster/bcl_cluster_distances_euclidean.h"
#include "model/bcl_model_retrieve_data_set_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> InputFeatures::s_Instance
    (
      GetObjectInstances().AddInstance( new InputFeatures())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    InputFeatures::InputFeatures() :
      m_DatasetRetriever(),
      m_RetrievedDataset(),
      m_Scale( false)
    {
    }

    //! @brief constructor taking members
    //! @param RETRIEVER method for retrieving feature vectors
    //! @param SCALE bool indicating if feature vectors should be scaled or not
    InputFeatures::InputFeatures( const util::Implementation< model::RetrieveDataSetBase> &RETRIEVER, bool SCALE) :
      m_DatasetRetriever( RETRIEVER),
      m_RetrievedDataset(),
      m_Scale( SCALE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new InputFeatures
    InputFeatures *InputFeatures::Clone() const
    {
      return new InputFeatures( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &InputFeatures::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief HandleInput gets data from input; puts it into the data construct with pointers to the objects
    //         in m_Objects
    //! @param IFSTREAM is the stream from which the input will be read
    //! @return returns data construct holding the distances between all objects for use in a LinkageInterface
    util::ShPtr
    <
      math::FunctionInterface
      <
        storage::VectorND< 2, util::SiPtr< const linal::Vector< float> > >, float
      >
    >
    InputFeatures::HandleInput( io::IFStream &IFSTREAM)
    {
      // retrieve the data set of features and results and set the member variable
      m_RetrievedDataset = m_DatasetRetriever->GenerateDataSet();

      // true if the data set needs to be scaled
      if( m_Scale)
      {
        // scale
        m_RetrievedDataset->GetFeaturesReference().CopyValues( Rescale( m_RetrievedDataset->GetFeaturesReference()));
      }

      linal::MatrixConstReference< float> features( m_RetrievedDataset->GetFeaturesReference());

      // make list that will hold the objects to be clustered
      util::ShPtr< storage::List< linal::Vector< float> > > cluster_objects( new storage::List< linal::Vector< float> >());

      // get the size of each feature
      const size_t feature_size( features.GetNumberCols());

      for
      (
        size_t feature_number( 0), number_features( features.GetNumberRows());
        feature_number < number_features;
        ++feature_number
      )
      {
        cluster_objects->PushBack( linal::Vector< float>( feature_size, features[ feature_number]));
      }

      SetInputObjects( cluster_objects);

      return util::ShPtr
        <
          math::FunctionInterface
          <
            storage::VectorND< 2, util::SiPtr< const linal::Vector< float> > >, float
          >
        >( new DistancesEuclidean< float>());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &InputFeatures::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DatasetRetriever, ISTREAM);
      io::Serialize::Read( m_RetrievedDataset, ISTREAM);
      io::Serialize::Read( m_Scale, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &InputFeatures::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DatasetRetriever, OSTREAM, INDENT);
      io::Serialize::Write( m_RetrievedDataset, OSTREAM, INDENT);
      io::Serialize::Write( m_Scale, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief rescales a matrix
    //! @param DATA the matrix to rescale
    //! @return matrix that has been rescaled
    linal::Matrix< float> InputFeatures::Rescale( const linal::MatrixConstInterface< float> &DATA)
    {
      // store rescale information here
      storage::Vector< math::Range< float> > rescale_ranges;
      const size_t feature_size( DATA.GetNumberCols());
      rescale_ranges.AllocateMemory( feature_size);

      // scale all features and results between 0 and 1, if possible
      math::Range< float> autoscale_range( 0.0, 1.0);

      // creating vector of ranges
      float tmp( 0);
      for( size_t col( 0); col < feature_size; ++col)
      {
        // set initial min and max
        float min( std::numeric_limits< float>::infinity());
        float max( -std::numeric_limits< float>::infinity());

        // iterate through all values in that feature col
        for( size_t row( 0), num_rows( DATA.GetNumberRows()); row < num_rows; ++row)
        {
          // set as min or max if qualifies
          tmp = DATA( row, col);
          if( tmp < min)
          {
            min = tmp;
          }
          else if( tmp > max)
          {
            max = tmp;
          }
        }
        // handle the case where min == max
        if( min == max)
        {
          min -= 0.5;
          max += 0.5;
        }
        // add to range vector
        rescale_ranges.PushBack( math::Range< float>( min, max));
      }

      // copy data so it can be rescaled
      linal::Matrix< float> rescaled_matrix( DATA);

      for( size_t row( 0), num_rows( DATA.GetNumberRows()); row < num_rows; ++row)
      {
        float *feature_row( rescaled_matrix[ row]);
        // rescale each feature
        for( size_t col( 0); col < feature_size; ++col)
        {
          feature_row[ col] = autoscale_range( DATA( row, col), rescale_ranges( col));
        }
      }

      return rescaled_matrix;
    }

  } // namespace cluster

} // namespace bcl
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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "cluster/bcl_cluster_input_pairwise_list.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
  } // namespace cluster
} // namespace bcl
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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "cluster/bcl_cluster_input_table.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
  } // namespace cluster
} // namespace bcl
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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "cluster/bcl_cluster_linkage_classes.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {

  } // namespace cluster

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< cluster::LinkageInterface< std::string, float> >, cluster::LinkageClasses< std::string, float> >;

  } // namespace util
} // namespace bcl
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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include forward header of this class
#include "cluster/bcl_cluster_node_description_from_file.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> NodeDescriptionFromFile::s_Instance
    (
      GetObjectInstances().AddInstance( new NodeDescriptionFromFile())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    NodeDescriptionFromFile::NodeDescriptionFromFile() :
      m_Descriptions()
    {
    }

    //! @brief constructor taking a filename and filling "m_Descriptions" based on the file's content
    //! @param FILENAME the path and name of the file which contains the strings and their descriptions
    NodeDescriptionFromFile::NodeDescriptionFromFile( const std::string &FILENAME) :
      m_Descriptions( GetDescriptions( FILENAME))
    {
    }

    //! @brief copy constructor
    //! @param DESCRIPTION the object which will be copied into this instance
    NodeDescriptionFromFile::NodeDescriptionFromFile
    (
      const NodeDescriptionFromFile &DESCRIPTION
    ) :
      m_Descriptions( DESCRIPTION.m_Descriptions)
    {
    }

    //! @brief virtual copy constructor
    //! @return new copy of this class
    NodeDescriptionFromFile *NodeDescriptionFromFile::Clone() const
    {
      return new NodeDescriptionFromFile( *this);
    }

    //! @brief virtual destructor
    NodeDescriptionFromFile::~NodeDescriptionFromFile()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &NodeDescriptionFromFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking a Wrapper<string> and returns a double which is the desription of the object
    //! @param MEMBER the object which was clustered for which a description is desired
    //! @return returns a double which is the desription of "MEMBER"
    double NodeDescriptionFromFile::operator()( const std::string &MEMBER) const
    {
      // create const iterator "itr" to the "MEMBER" key in "m_Descriptions"
      storage::Map< std::string, double>::const_iterator itr( m_Descriptions.Find( MEMBER));

      // make sure that "itr" is not at the end of "m_Descriptions" - "MEMBER" could not be found in "m_Descriptions"
      if( itr == m_Descriptions.End())
      {
        // write a warning message that "MEMBER" could not be found
        BCL_MessageStd( MEMBER + " could not be found in descriptor list");

        // return undefined double
        return util::GetUndefined< double>();
      }

      // return the value description that corresponds to "MEMBER"
      return itr->second;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &NodeDescriptionFromFile::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Descriptions, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &NodeDescriptionFromFile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Descriptions, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief GetDescriptions takes a file and reads in the descriptions for the clustered objects
    //! @param "FILENAME" the path and name of the file which contains the strings and their descriptions
    //! @return return a map which contains the clustered string objects and each numberical description
    storage::Map< std::string, double> NodeDescriptionFromFile::GetDescriptions
    (
      const std::string &FILENAME
    ) const
    {
      // create IFStream "read"
      io::IFStream read;

      // open "read" and bind it to "FILENAME"
      io::File::MustOpenIFStream( read, FILENAME);

      // create Vector of Vectors of strings "split_lines" and initialize with the split lines in "FILENAME"
      const storage::Vector< storage::Vector< std::string> > split_lines
      (
        util::SplittedStringLineListFromIStream( read)
      );

      // create Map "descriptions" which will store all of the descriptions and clustered strings in "FILENAME"
      storage::Map< std::string, double> descriptions;

      // iterate through "split_lines" in order to fill "descriptions"
      for
      (
        storage::Vector< storage::Vector< std::string> >::const_iterator
          itr( split_lines.Begin()), itr_end( split_lines.End());
        itr != itr_end;
        ++itr
      )
      {
        // create const size_t "line_num_cols" and initialize with the number of columns in
        // the current line denoted by "itr"
        const size_t line_num_cols( itr->GetSize());

        // true if there are less than two columns in the current line of "FILENAME" denoted by "itr"
        if( line_num_cols < 2)
        {
          // write a warning message that there is not at least 2 columns in current line
          BCL_MessageStd
          (
            "Only " + util::Format()( line_num_cols) + " on line " +
            util::Format()( itr - split_lines.Begin()) + " of file " + FILENAME
          );

          // go to next line
          continue;
        }

        // true if there are more than two columns in the current line of "FILENAME" denoted by "itr"
        // any addition columns will be ignored
        else if( line_num_cols > 2)
        {
          // write a warning message that additional columns will be ignored
          BCL_MessageStd
          (
            "More than 2 columns. Additional columns in " + FILENAME + " will be ignored"
          );
        };

        // create const reference to string "string_object" and initialize with the string object denoted by "itr"
        const std::string &string_object( itr->operator()( 0));

        // create const double "string_object" and initialize with the double description denoted by "itr"
        const double double_description( util::ConvertStringToNumericalValue< double>( itr->operator()( 1)));

        // create pair "insertion" and initialize with the results of trying to
        // insert the string and its numerical description into "descriptions" from the first and second columns
        std::pair< storage::Map< std::string, double>::const_iterator, bool> insertion
        (
          descriptions.Insert
          (
            std::pair< std::string, double>
            (
              string_object, double_description
            )
          )
        );

        // true if the insertion of the string and its description was unsuccessful
        if( !insertion.second)
        {
          // write warning message that the insertion failed maybe because that key already exists in "descriptions"
          BCL_MessageStd
          (
            "Failed to insert " + string_object + " and its description " +
            util::Format()( double_description) + " Maybe it already has a description"
          );
        }
      }

      // return "descriptions"
      return descriptions;
    }

  } // namespace cluster
} // namespace bcl
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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "cluster/bcl_cluster_output_classes.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {

  } // namespace cluster

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< cluster::OutputInterface< std::string, float> >, cluster::OutputClasses< std::string, float> >;

  } // namespace util
} // namespace bcl
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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "cluster/bcl_cluster_output_pymol_label_string.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
  } // namespace cluster
} // namespace bcl
