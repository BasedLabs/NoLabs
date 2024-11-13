export type AminoAcid = {
    name: string;
    sequence: string;
    go: {
        [name: string]: {
            name: string,
            namespace: string,
            edges: { [name: string]: Array<string> }
        }
    };
};

