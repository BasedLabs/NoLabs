export {}

declare module '@vue/runtime-core' {
  interface ComponentCustomProperties {
    $loaderShow: (message: string) => void;
  }
}
