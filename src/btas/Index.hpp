#ifndef __BTAS_INDEX_HPP
#define __BTAS_INDEX_HPP 1

//
//  [1] I want to write tensor contraction as
//       + contract(1.0, A, "al,px,rx", B, "rx,qx,br", 1.0, C, "al,px,qx,br");
//
//      Function is declared as
//       + void contract(double, Tensor<double>&, string, Tensor<double>&, string, Tensor<double>&, string);
//
//  [2] I also want to write QN merging as
//       + merge(A, "i,j,k,l", B, "[i,j],k,l");
//       + merge(A, "i,j,k,l", B, "k,l,m,n", AM, "i,j,[k,l]", BM, "[k,l],m,n");
//
//  Issues:
//   + What objects are involved in "Index" class?
//

namespace btas {

class Label {

public:

   Label() { }

   virtual ~Label() { }

   void reset(const std::string& expr)
   {
      label_ = expr;
   }

   template<typename _T>
   void reset(const _T& x)
   {
      std::ostringstream expr; expr << x;
      label_ = expr.str();
   }

   const std::string& str() const { return label_; }

   virtual bool is_group() const { return false; }

private:

   std::string label_;

};



class Index : public Label {

public:

   Index() { }

  ~Index() { }

   Index(const string& expr) { reset(expr); }

   template<typename... Args>
   void reset(const Args&... args)
   {
      index_.clear();
      __reset_from_args(args);
   }

   void reset(const string& expr)
   {
      index_.clear();
      __reset_from_string(expr);
   }

   bool is_group() const { return !index_.empty(); }

   const std::string& operator[] (size_t n) const { index_[n]->str(); }

private:

   void __reset_from_string(const std::string& expr)
   {
      std::vector<std::string> labs = parse(expr);
      for(size_t i = 0; i < labs.size(); ++i) {
         if(labs[i].find(',') != std::string::npos)
            index_.push_back(std::unique_ptr<Label>(new Index(labs[i])));
         else
            index_.push_back(std::unique_ptr<Label>(new Label(labs[i])));
      }
   }

   template<typename Arg1, typename... Args>
   void __reset_from_args(const Arg1& expr, const Args&... rest)
   {
      index_.push_back(std::unique_ptr<Label>(new Label(expr)));
      __reset_from_args(rest);
   }

   template<typename... Args>
   void __reset_from_args(const std::string& expr, const Args&... rest)
   {
      if(expr.find(',') != std::string::npos)
         index_.push_back(std::unique_ptr<Label>(new Index(expr)));
      else
         index_.push_back(std::unique_ptr<Label>(new Label(expr)));
      __reset_from_args(rest);
   }

   void __reset_from_args() { }

   std::vector<std::unique_ptr<Label>> index_;

};

std::vector<std::string> parse(const std::string& expr)
{
   std::vector<std::string> labs;
   for(size_t i = 0; i < expr.size(); ++i) {
      if(expr[i] == '[') {
         size_t j = expr.find(']', i);
         BTAS_RUNTIME_ASSERT(j != std::string::npos, "found '[' with no closure ']'");
         size_t k = expr.find('[', i);
         BTAS_RUNTIME_ASSERT(k == std::string::npos || k > j, "found another '[' before closure ']'");
         labs.push_back(std::string(expr, i+1, j-i-1));
         i = j;
      }
      else {
         size_t j = expr.find(',', i);
         if(j == std::string::npos) j = expr.size();
         labs.push_back(std::string(expr, i, j-i-1));
         i = j;
      }
   }
   return labs;
}

std::vector<Index> make_index(const std::string& expr)
{
   std::vector<std::string> labs = parse(expr);
   std::vector<Index> indxs;
   indxs.reserve(labs.size());
   for(size_t i = 0; i < labs.size(); ++i) {
      indxs.push_back(Index(labs[i]));
   }
   return indxs;
}

} // namespace btas

#endif // __BTAS_INDEX_HPP
