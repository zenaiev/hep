void hello() {
  printf("Hello world\n");
  TH1F hist("hist", "My histo", 300, 0., 3.);
  hist.Fill(1.2);
  hist.Print();
}