#ifndef STATISMO_ASMIMAGEPREPROCESSOR_H
#define STATISMO_ASMIMAGEPREPROCESSOR_H

#include "CommonTypes.h"

namespace statismo {

    template <typename TPointSet>
    class ASMPreprocessedImage {
    public:
        typedef typename Representer<TPointSet>::PointType PointType;
        virtual ~ASMPreprocessedImage() {}
        virtual bool IsDefined(const PointType& point) const = 0;
        virtual VectorType Evaluate(const PointType& point) const = 0;
    };

    template <typename TPointSet, typename TImage>
    class ASMImagePreprocessor {
    public:
        virtual ~ASMImagePreprocessor() {}
        virtual ASMPreprocessedImage<TPointSet>* Preprocess(const TImage* image) const = 0;
    };


    template<typename TPointSet, typename TImage>
    class ASMImagePreprocessorFactory {
        typedef ASMImagePreprocessorFactory<TPointSet, TImage> ASMImagePreprocessorFactoryType;
    private:
        static std::vector<const ASMImagePreprocessorFactoryType*> *implementations() {
            static std::vector<const ASMImagePreprocessorFactoryType*> impls;
            return &impls;
        }
        ASMImagePreprocessorFactory(const ASMImagePreprocessorFactoryType& o) { }
        ASMImagePreprocessorFactory& operator=(const ASMImagePreprocessorFactoryType& o) {}

    protected:
        ASMImagePreprocessorFactory() {}

    public:
        virtual std::string GetDescriptor() const = 0;
        virtual const ASMImagePreprocessor<TPointSet, TImage>* Instantiate(const H5::Group& h5Group) const = 0;

        static std::vector<const ASMImagePreprocessorFactoryType*> GetImplementations() {
            return *implementations();
        }

        static void RegisterImplementation(const ASMImagePreprocessorFactoryType* impl) {
            if (GetImplementation(impl->GetDescriptor())) {
                //ignoring already-registered implementation
                return;
            }
            implementations()->push_back(impl);
        }

        static const ASMImagePreprocessorFactoryType* GetImplementation(std::string descriptor) {
            std::vector<const ASMImagePreprocessorFactoryType* > impls = GetImplementations();
            for (typename std::vector<const ASMImagePreprocessorFactoryType* >::iterator impl = impls.begin(); impl != impls.end(); ++impl) {
                if ((*impl)->GetDescriptor() == descriptor) {
                    return *impl;
                }
            }
            return 0;
        }
    };

}
#endif //STATISMO_ASMIMAGEPREPROCESSOR_H
