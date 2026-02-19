use proc_macro::TokenStream;
use quote::quote;
use syn::parse_macro_input;
use syn::ItemFn;

#[proc_macro_attribute]
pub fn profiling(_attr: TokenStream, item: TokenStream) -> TokenStream {
    let input = parse_macro_input!(item as ItemFn);

    // If the "profiling" feature isn't enabled in the macro crate,
    // just spit the function back out exactly as it was.
    if !cfg!(feature = "profiling") {
        return quote!(#input).into();
    }

    let vis = &input.vis;
    let sig = &input.sig;
    let block = &input.block;
    let name = sig.ident.to_string();

    quote! {
        #vis #sig {
            crate::profiling::push(#name);
            let result = (move || #block)();
            crate::profiling::pop();
            result
        }
    }
    .into()
}
