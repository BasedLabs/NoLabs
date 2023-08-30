describe('Main page', () => {
  const inputAminoAcidSequenceId = '#inputSequence';
  const submitAminoAcidSequenceId = '#submitInference';
  const localisationImageClass = '.localisation-image';
  const geneOntologySvg = '#geneOntologyContainer svg';

  it('Opened main page', () => {
    cy.visit('http://127.0.0.1:5000');
  });

  it('Amino-acid sequence input is visible', () => {
    cy.visit('http://127.0.0.1:5000');
    cy.get(inputAminoAcidSequenceId).should('exist');
  });

  it('Inference result is visible', () => {
    cy.visit('http://127.0.0.1:5000');
    cy.get(inputAminoAcidSequenceId).click();
    cy.get(inputAminoAcidSequenceId).type(`MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA`);
    cy.get(submitAminoAcidSequenceId).click();
    cy.get(localisationImageClass).should('exist');
    cy.get(geneOntologySvg, { timeout: 10000 }).should('exist');
  });
})